from botorch.sampling import SobolQMCNormalSampler
from botorch.acquisition.objective import GenericMCObjective
from botorch.models import FixedNoiseGP, SingleTaskGP
from botorch.models.multitask import KroneckerMultiTaskGP
from botorch.models.transforms.outcome import Standardize
from botorch.fit import fit_gpytorch_model
from botorch.optim import optimize_acqf
from botorch.optim.fit import fit_gpytorch_mll_torch
from botorch.acquisition import qExpectedImprovement
from botorch.generation import MaxPosteriorSampling
from botorch.utils.transforms import normalize, unnormalize

from gpytorch.mlls import ExactMarginalLogLikelihood
from gpytorch.constraints import GreaterThan
from gpytorch.likelihoods import GaussianLikelihood, MultitaskGaussianLikelihood

from botorch.test_functions.synthetic import Rosenbrock, Ackley

import torch
from torch.quasirandom import SobolEngine
from torch.optim.adam import Adam
from functools import partial

import subprocess
import numpy as np

import json
import math

dtype = torch.double
device = torch.device("cpu")
    
def fit_SingleTaskGP(X, Y):
    model = SingleTaskGP(
        X, Y,
        likelihood=GaussianLikelihood(noise_constraint=GreaterThan(1e-6)),
        outcome_transform=Standardize(m=Y.shape[-1]) # Standardize Y for training, then untransform Y for predictions
    )
    mll = ExactMarginalLogLikelihood(model.likelihood, model)
    fit_gpytorch_model(mll)
    
    return model

def fit_KroneckerMultiTaskGP(X, Y):
    model = KroneckerMultiTaskGP(
        X, Y,
        likelihood=MultitaskGaussianLikelihood(num_tasks=Y.shape[-1], noise_constraint=GreaterThan(1e-6)),
        outcome_transform=Standardize(m=Y.shape[-1]) # Standardize Y for training, then untransform Y for predictions
    )

    mll = ExactMarginalLogLikelihood(model.likelihood, model)
        
    fit_gpytorch_mll_torch(
        mll,
        optimizer=partial(Adam, lr=0.1),
        step_limit=3000,
    )
    
    return model

def get_turbo_bounds(X, Y, length):
    # Scale the TR to be proportional to the lengthscales
    x_center = X[objective(Y).argmax(), :].clone()

    # hack for single and multi task GPs
    try:
        weights = model.covar_module.base_kernel.lengthscale.squeeze().detach()
    except Exception as e:
        try:
            weights = model.covar_module.data_covar_module.lengthscale.squeeze().detach()
        except Exception as e:
            pass
    
    # hack for 1-D problems
    if len(weights.shape) == 0: weights = weights.unsqueeze(-1)
                    
    weights = weights / weights.mean()
    weights = weights / torch.prod(weights.pow(1.0 / len(weights)))
    tr_lb = torch.clamp(x_center - weights * length / 2.0, 0.0, 1.0)
    tr_ub = torch.clamp(x_center + weights * length / 2.0, 0.0, 1.0)

    return tr_lb, tr_ub


def generate_batch_ExpectedImprovement(X, Y, model, objective, batch_size, tr_lb, tr_ub):
    sampler = SobolQMCNormalSampler(1024, collapse_batch_dims=True)

    ei = qExpectedImprovement(model, best_f=objective(Y).max(), objective=objective, sampler=sampler, maximize=True)
    
    X_batch, acq_value = optimize_acqf(
        ei,
        bounds=torch.stack([tr_lb, tr_ub]),
        q=batch_size,
        num_restarts=10,
        raw_samples=512
    )
    
    return X_batch

def generate_batch_ThompsonSampling(X, Y, model, objective, batch_size, tr_lb, tr_ub):
    n_candidates = 10000

    dim = X.shape[-1]

    sobol = SobolEngine(dim, scramble=True)
    pert = sobol.draw(n_candidates).to(dtype=dtype, device=device)
    pert = tr_lb + (tr_ub - tr_lb) * pert

    # Create a perturbation mask
    prob_perturb = min(20.0 / dim, 1.0)
    mask = torch.rand(n_candidates, dim, dtype=dtype, device=device) <= prob_perturb
    ind = torch.where(mask.sum(dim=1) == 0)[0]

    if len(ind) > 0:
        mask[ind, torch.randint(0, dim - 1, size=(len(ind),), device=device)] = 1
        
    # Create candidate points from the perturbations and the mask
    x_center = X[objective(Y).argmax(), :].clone()
    X_cand = x_center.expand(n_candidates, dim).clone()
    X_cand[mask] = pert[mask]

    # Sample on the candidate points
    thompson_sampling = MaxPosteriorSampling(model=model, replacement=False, objective=objective)
    
    with torch.no_grad():  # We don't need gradients when using TS
        X_batch = thompson_sampling(X_cand, num_samples=batch_size)

    return X_batch

def load():
    try:
        with open('state.json', 'r') as fd:
            return json.load(fd)
    except FileNotFoundError:
        return None

def save(dict):
    with open('state.json', 'w') as fd:
        json.dump(dict, fd)

# =============================================================================
# MAIN
# =============================================================================
torch.cuda.empty_cache()

dim = 4

# R function
bounds = torch.tensor([[0.0 for i in range(dim)], [1.0 for i in range(dim)]], dtype=dtype, device=device)

def exec(command):
    return subprocess.Popen(command, shell = True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True, text=True) #, capture_output = capture)#command) #, shell = True, capture_output = True)

def evaluate(X_batch):
    np.savetxt('X.txt', X_batch.cpu().numpy())
    p = subprocess.Popen("bash ./run_R.sh", shell = True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True, text=True)
    outs, errs = p.communicate(timeout=None)
    print(outs)
    print(errs)
    p.wait()
    
    Y_batch = np.loadtxt('Y.txt', ndmin=2)
    Xout = []
    Yout = []
    for i in range(0, len(X_batch)):
        for j in range(len(Y_batch[i])):
            if math.isnan(Y_batch[i][j]):
                continue
            Xout.append(X_batch[i].tolist())
            Yout.append(Y_batch[i].tolist())
    return torch.tensor(list(Xout), dtype=dtype, device=device), torch.tensor(list(Yout), dtype=dtype, device=device)

init_samples = 100 # the more parameters the more the samples 500
init_batches = 1
max_evaluations = 10000
batch_size = 64 # how many samples for every iteration
objective = GenericMCObjective(lambda yi, X: torch.sum(yi, dim=-1)) # function to convert a multi-objective output to a single value; a simple sum by default

# Turbo parameters
failure_counter = 0
failure_tolerance = 3 # setting this to torch.inf will never shrink the search range
turbo_length = 1.0 # start with 1.0 to cover the whole unit-hypercube and progressively shrink it

# Best value so far
Y_max: float = -torch.inf

# Try to load from previous state
state = load()
if state is not None:
    turbo_length, Y_max, X, Y = state['turbo_length'], state['Y_max'], torch.tensor(state['X'], dtype=dtype, device=device), torch.tensor(state['Y'], dtype=dtype, device=device)
else:
    # Evaluate initial samples
    X = torch.tensor([], dtype=dtype, device=device)
    Y = torch.tensor([], dtype=dtype, device=device)
    
    sobol = SobolEngine(dimension=dim, scramble=True)
    X_init = sobol.draw(n=init_samples).to(dtype=dtype, device=device)
    
    # Evaluate init samples in batches to avoid large array jobs in Slurm
    init_batch_size = int(len(X_init)/init_batches) 
    for batch in range(0, init_batches):
        batch_start = batch * init_batch_size
        batch_end = batch_start + init_batch_size
    
        if batch == init_batches-1:
            batch_end = batch_start + len(X_init) - (init_batches-1) * init_batch_size
    
        X_batch, Y_batch = evaluate(X_init[batch_start:batch_end])
    
        # Add init samples to the dataset
        X = torch.cat((X, X_batch.to(dtype=dtype, device=device)), axis=0)
        Y = torch.cat((Y, Y_batch.to(dtype=dtype, device=device)), axis=0)
    
    save({'turbo_length': turbo_length, 'Y_max': Y_max, 'X': X.cpu().numpy().tolist(), 'Y': Y.cpu().numpy().tolist()})

print(f'{len(X)}) Turbo: {turbo_length}, Y_max = {objective(Y).max().item()} at {unnormalize(X[objective(Y).argmax()], bounds).squeeze().cpu().numpy()}')

while len(X) < max_evaluations:
    # Fit model
    if Y.shape[-1] == 1:
        model = fit_SingleTaskGP(X, Y) # single objective output
    else:
        model = fit_KroneckerMultiTaskGP(X, Y) # multi-objective output

    # Increase failure_counter if no improvement to the lossfunction
    if objective(Y).max().item() > Y_max: 
        failure_counter = 0
        Y_max = objective(Y).max().item()
    else:
        failure_counter += 1

    # Shrink Turbo search range when failure_counter = failure_tolerance
    if failure_counter == failure_tolerance:
        turbo_length /= 2.0
        failure_counter = 0

    # Calculate Turbo search range using the covariance matrix
    tr_lb, tr_ub = get_turbo_bounds(X, Y, turbo_length)

    # Generate next batch of points
    X_batch_ei = generate_batch_ExpectedImprovement(X, Y, model, objective, batch_size, tr_lb, tr_ub)
    #X_batch_ts = generate_batch_ThompsonSampling(X, Y, model, objective, batch_size, tr_lb, tr_ub)
    # X_batch = torch.cat((X_batch_ei, X_batch_ts), dim=0)

    # Evaluate batch
    X_batch, Y_batch = evaluate(X_batch_ei)
    
    # Append new evaluations to dataset
    X = torch.cat((X, X_batch), dim=0)
    Y = torch.cat((Y, Y_batch), dim=0)

    print(f'{len(X)}) Turbo: {turbo_length}, Y_max = {Y_max} at {unnormalize(X[objective(Y).argmax()], bounds).squeeze().cpu().numpy()}')

    save({'turbo_length': turbo_length, 'Y_max': Y_max, 'X': X.cpu().numpy().tolist(), 'Y': Y.cpu().numpy().tolist()})
