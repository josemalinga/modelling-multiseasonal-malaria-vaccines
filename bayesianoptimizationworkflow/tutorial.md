---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.15.2
  kernelspec:
    display_name: Python [conda env:miniconda3]
    language: python
    name: conda-env-miniconda3-py
---

# Botorch Quick Tutorial


This is a quick and simple Botorch tutorial. This will cover how to:
* Manipulate Pytorch tensors
* Fit a gaussian process
* Use an acquisition function
* Do Bayesian optimization with adaptive sampling

See: https://botorch.org/docs/introduction


## Ackley

We start by importing a few modules:
* Botorch provides several test functions that we can use as examples. Let us use Ackley.
* We also import the unnormalize function, which does exactly what it means
* matplotlib is a standard plotting library in Python
* torch is part of Pytorch and is used to manipulate tensors

```python
from botorch.test_functions.synthetic import Ackley
from botorch.utils.transforms import unnormalize
import matplotlib.pyplot as plt # To create plots
import torch # to manipulate Pytorch tensors
```

Now we are going to generate some points (x, y) using the Ackley function provided by Botorch. We will only use `dim=1` for now. Note that the function expects a list of list (or rather a 2D tensor) of the form:
```
[
    [x0_0, x0_1, ..., x0_m],
    [x1_0, x1_1, ..., x1_m],
    ...
    [xn_0, xn_1, ..., xn_m]
]
```
where `m` is the number of dimensions. However since we have only one dimension for now we need to pass:
```
[
    [x1],
    [x2],
    ...
    [xn] 
]
```
We generate `100` points between `0` and `1` and use the `unsqueeze(-1)` method to add one inner dimension to the data to conform to the expectation of Botorch. Note that we reverse the operation by using the `squeeze(-1)` method. We then sample the Ackley function between `-4` and `4` thanks to the unnormalize function. Here it is also a good time to `unsqueeze` `y` to match the shape of `x`.

**Note 1: It is good practice to work with data between 0 and 1 and only unnormalize when we need to evaluate the function.**

**Note 2: It is common to have to adjust the shape of the data using `squeeze()` and `unsqueeze()` while using Botorch.**

```python
ackley_function = Ackley(dim=1)

bounds = [-4, 4]
X = torch.linspace(0, 1, 100).unsqueeze(-1)
Y = ackley_function(unnormalize(X, bounds)).unsqueeze(-1)

#print(x)
#print(y)
plt.plot(unnormalize(X, bounds), Y)
```

## Fit a Gaussian process

Now that we have our `x` and `y` data, we can fit a gaussian process. Botorch provides multiple model implementations, see: https://botorch.org/docs/models.

We will use the simple `SingleTaskGP` model with the following parameters:
* `likelihood=GaussianLikelihood(noise_constraint=GreaterThan(1e-4))` is the standard likelihood for regression with support for missing values. Assumes a standard homoskedastic noise model (inferred noise). We can provide a constrain on the noise, but it can be finnicky with low or high minimum / maxium noise constraints.
*  `outcome_transform=Standardize(m=y.shape[-1])` will automatically standardize the `y` data before training the model, and unstdandardize it when using the model to make predictions.

**Note: model fitting works better with standardized data, so it's good practice to add this paramater. Even better here the whole process is transparent to us.**

```python
from botorch.models import SingleTaskGP
from botorch.models.transforms.outcome import Standardize

from gpytorch.constraints import GreaterThan
from gpytorch.likelihoods import GaussianLikelihood

model = SingleTaskGP(
    X, Y,
    likelihood=GaussianLikelihood(noise_constraint=GreaterThan(1e-4)),
    outcome_transform=Standardize(m=Y.shape[-1])
)
```

Now that we have created our model, we can fit it to our data. 
* Here "fitting" refers to the optimizing the hyperparameters of the underlying kernel function
* The kernel function is used to compute the underlying covariance matrix
* Here a Matern kernel (5/2) is used by default, which works well in most cases
* From the covariance matrix, we can then make predictions using the Gaussian Likelihood

In order to fit the hyperparameters of the Gaussian process, we extract and fit the marginal log-likelihood using a function provided by gpytorch. If you want to learn more about this whole process, you can refer to the official documentation of Botorch and Pytorch.

```python
from botorch.fit import fit_gpytorch_model

from gpytorch.mlls import ExactMarginalLogLikelihood

mll = ExactMarginalLogLikelihood(model.likelihood, model)
fit_gpytorch_model(mll)
```

## Use the Gaussian Process


Our model is now fitted and ready to be used. We can use it to make predictions of the mean and variance. We generate some points and we sample the posterior. Here we use `torch.no_grad()` to tell pytorch that we don't need to perform gradient calculation.

```python
X_model = torch.linspace(0, 1, 100).unsqueeze(-1) # generate 100 equally spaced points on [0, 1]

# sample the posterior
with torch.no_grad():
    posterior = model.posterior(x_model)

Y_model_mean = posterior.mean # predicted mean
Y_model_std = torch.sqrt(posterior.variance) # standard deviation from the predicted variance
Y_model_low = Y_model_mean - 1.96 * Y_model_std # 95% low
Y_model_high = Y_model_mean + 1.96 * Y_model_std # 95% high

plt.figure(figsize=(8,6))
plt.plot(unnormalize(X, bounds), Y, linewidth=6, label="Original function")
plt.plot(unnormalize(X_model, bounds), Y_model_mean, linewidth=3, linestyle='--', label="Predicted mean")
plt.fill_between(unnormalize(X_model, bounds).squeeze(), Y_model_low.squeeze(), Y_model_high.squeeze(), alpha=0.5, label="95% CI")
plt.legend()
```

Well, for such a simple one-dimensional function, the model is matching the original data perfectly!


# Bayesian optimization with adaptive sampling


Suppose that we don't know the shape of the Ackley function, but we would like to find its minimum value, and more importantly the `x` that minimizes the function. It could be any `x` between `[-4,4]`. Here we can argue that the process is straightforward because there is only one dimension and the function can be evaluated very quickly. But in practice, Bayesan optimization is useful when it takes a lot of time and resources to evaluate the function as it typically requires less sample to converge towards the optimal solution compared to other optimization methods.


## Init samples

In order to boostrap the whole process, we need at least a few samples. We are going to sample `2` points.

**Note: most approach maximize by default. Because in this example we want to find the minimum point, we negate the function and work with negative values instead.**

```python
from torch.quasirandom import SobolEngine

dim = 1
n_init = 2
sobol = SobolEngine(dimension=dim, scramble=True)

ackley_function = Ackley(dim=dim, negate=True)

bounds = [-5, 5]
X = sobol.draw(n=n_init) # the data is alreay in the correct shape, no need to unsqueeze
Y = ackley_function(unnormalize(X, bounds)).unsqueeze(-1)

#print(X)
#print(Y)
plt.scatter(unnormalize(X, bounds), Y)
```

## Fit Gaussian process

Not so easy to find the minimum now. But this is enough to boostrap the optimization process. We first define a function to fit our Gaussian process.

```python
def fit_SingleTaskGP(X, Y):
    model = SingleTaskGP(
        X, Y,
        likelihood=GaussianLikelihood(noise_constraint=GreaterThan(1e-6)),
        outcome_transform=Standardize(m=Y.shape[-1]) # Standardize Y for training, then untransform Y for predictions
    )
    mll = ExactMarginalLogLikelihood(model.likelihood, model)
    fit_gpytorch_model(mll)
    
    return model
```

And we fit our model to our two data points (which is not great, but it's a start).

```python
model = fit_SingleTaskGP(X, Y)
```

## Plot predictions

We now define a function to show the prediction of the Gaussian process. This will be useful to monitor the progress.

```python
def plot_predictions(X, Y, model):
    x_model = torch.linspace(0, 1, 100).unsqueeze(-1)
    with torch.no_grad():
        posterior = model.posterior(x_model)
    y_model_mean = posterior.mean
    y_model_std = torch.sqrt(posterior.variance)
    y_model_low = y_model_mean - 1.96 * y_model_std
    y_model_high = y_model_mean + 1.96 * y_model_std
    
    plt.figure(figsize=(8,6))
    plt.plot(unnormalize(x_model, bounds), y_model_mean, linewidth=3, linestyle='--', color='C1', label="Predicted mean")
    plt.fill_between(unnormalize(x_model, bounds).squeeze(), y_model_low.squeeze(), y_model_high.squeeze(), alpha=0.3, color='C1', label="95% CI")
    plt.scatter(unnormalize(X, bounds), Y, linewidth=6, marker='o', label="Observed data")
    plt.ylim([0, -20])
    plt.legend()
```

We can see that the fitted model is a ble to predict the known data points and provide some uncertainty bounds everywhere else.

```python
plot_predictions(X, Y, model)
```

## Acquisition function

Ok, as expected there is a lot of uncertainty at this point, we can use this to our advantage. The entire point of doing Bayesian optimization with Gaussian processes is to use the surrogate model to help us explore the space efficiently. This is an exploitation vs exploration dilemma, also known as multi-armed bandit problem: should we evaluate the function where the prediction is at its minimum? Or should we evaluate the function where the uncertainty is high with the hope of finding a better value?

There are multiple ways to solve this dilemma. In Bayesian optimization, the standard approach is to use an acquisition function such as `Expected Improvement`, which will return the best candidate points based on the predicted mean and variance.

Botorch provide such acquisition functions, let us define a function that will return a new batch of points using `ExpectedImprovement`.

We use the Monte-Carlo version of the acquisition function because it can return more than 1 point if needed.

```python
from botorch.sampling import SobolQMCNormalSampler
from botorch.acquisition import qExpectedImprovement
from botorch.optim import optimize_acqf

def generate_batch_ExpectedImprovement(X, Y, model, batch_size, tr_lb, tr_ub):
    sampler = SobolQMCNormalSampler(1024, collapse_batch_dims=True)

    ei = qExpectedImprovement(model, best_f=Y.max(), sampler=sampler)
    
    X_batch, acq_value = optimize_acqf(
        ei,
        torch.stack([tr_lb, tr_ub]),
        q=batch_size,
        num_restarts=10,
        raw_samples=512
    )
    
    return X_batch
```

The function returns the best candidate point given the mean and variance of the surrogate model. Note that there are other parameters such as the number of samples and number of restarts that can be increased in order to optimize the acquisition function more thoroughly at the cost of additional time and memory.

```python
batch_size = 1
X_batch = generate_batch_ExpectedImprovement(X, Y, model, batch_size, torch.tensor([0.]), torch.tensor([1.]))
unnormalize(X_batch, bounds)
```

We can now evaluate the candidate(s) and plot them.

```python
Y_batch = ackley_function(unnormalize(X_batch, bounds)).unsqueeze(-1)

#print(X_batch)
#print(Y_batch)
plot_predictions(X, Y, model)
plt.scatter(unnormalize(X_batch, bounds), Y_batch, color='red', label="New point")
plt.legend()
```

At this point, ExpectedImprovement will often return a point where the uncertainty is very high. After evaluating the point, we will gain new knowledge on the function! We can then re-fit the Gaussian process by including the new data-point and observe the new predictions. Uncertainty is lower around the new point. Perhaps we have not yet found a better value, but we have improved the surrogate model.

```python
# concatenate the new batch points to our X and Y data
X = torch.cat((X, X_batch), dim=0)
Y = torch.cat((Y, Y_batch), dim=0)

model = fit_SingleTaskGP(X, Y)
plot_predictions(X, Y, model)
plt.scatter(unnormalize(X_batch, bounds), Y_batch, color='red', label="New point")
plt.legend()
```

## Optimization loop

All we have to do now is repeat the last steps in a loop until we converge on an optimal value.

```python
# Init samples
dim = 1
n_init = 2
batch_size = 2

sobol = SobolEngine(dimension=dim, scramble=True)

ackley_function = Ackley(dim=dim, negate=True)

bounds = [-4, 4]
X = sobol.draw(n=n_init) # the data is alreay in the correct shape, no need to unsqueeze
Y = ackley_function(unnormalize(X, bounds)).unsqueeze(-1)
```

```python
# Loop (ctrl+enter to run the same cell)
model = fit_SingleTaskGP(X, Y)

X_batch = generate_batch_ExpectedImprovement(X, Y, model, batch_size, torch.tensor([0.]), torch.tensor([1.]))
Y_batch = ackley_function(unnormalize(X_batch, bounds)).unsqueeze(-1)

plot_predictions(X, Y, model)
plt.scatter(unnormalize(X_batch, bounds), Y_batch, color='red', label="New point")
plt.legend()

X = torch.cat((X, X_batch), dim=0)
Y = torch.cat((Y, Y_batch), dim=0)

print(f'Evaluations: {len(Y)}, Minimum: {-Y.max()}')
```

## Multiple dimensions

We can change the example above to use an arbitrary number of input dimensions. We define a new function to plot the predictions in one dimension only, along a cut in the middle of the unit-hypercube (at x = [0.5, 0.5, ..., 0.5] except for one variable dimension. We also include the observed data that fall within the size of the cut (this is not great but it is still useful to observe progress).

```python
def plot_predictions_cut(X, Y, model, at=0.5, size=0):
    fig = plt.figure(figsize=(30,4))

    for d in range(X.shape[-1]):
        fig.add_subplot(1,X.shape[-1],d+1)
        
        cut = torch.linspace(0, 1, 100)
        X_model = torch.zeros(size=(100,X.shape[-1])) + at
        X_model[:,d] = cut
    
        X_inside_cut = []
        Y_inside_cut = []
        for i in range(len(X)):
            add = True
            for j in range(X.shape[-1]):
                if j != d:
                    if X[i,j] < at-size or X[i,j] > at+size:
                        add = False
                        break
            if add:
                X_inside_cut.append(X[i].tolist())
                Y_inside_cut.append(Y[i].tolist())
        X_inside_cut = torch.tensor(X_inside_cut)
        Y_inside_cut = torch.tensor(Y_inside_cut)
    
        with torch.no_grad():
            posterior = model.posterior(X_model)
            
            Y_model_mean = posterior.mean
            Y_model_std = torch.sqrt(posterior.variance)
            Y_model_low = Y_model_mean - 1.96 * Y_model_std
            Y_model_high = Y_model_mean + 1.96 * Y_model_std
    
            plt.plot(unnormalize(X_model, bounds).squeeze()[:,d], Y_model_mean.squeeze(), color='C1', label="Predicted mean")
            plt.fill_between(unnormalize(X_model, bounds).squeeze()[:,d], Y_model_low.squeeze(), Y_model_high.squeeze(), alpha=0.3, color='C1', label="95% CI")
            if len(X_inside_cut) > 0:
                plt.scatter(unnormalize(X_inside_cut, bounds)[:,d].squeeze(), Y_inside_cut.squeeze(), marker='o', label="Observed data")
            # plt.ylim([0, -20])
            plt.title(f'd={d}')
    plt.legend()
```

```python
# Init samples
dim = 5
n_init = 100
batch_size = 5

sobol = SobolEngine(dimension=dim, scramble=True)

bounds = [-4, 4]
ackley_function = Ackley(dim=dim, negate=True)

X = sobol.draw(n=n_init) # the data is alreay in the correct shape, no need to unsqueeze
Y = ackley_function(unnormalize(X, bounds)).unsqueeze(-1)
```

```python
# Loop (ctrl+enter to run the same cell)
model = fit_SingleTaskGP(X, Y)

tr_lb = torch.zeros(X.shape[-1])
tr_ub = torch.ones(X.shape[-1])
X_batch = generate_batch_ExpectedImprovement(X, Y, model, batch_size, tr_lb, tr_ub)
Y_batch = ackley_function(unnormalize(X_batch, bounds)).unsqueeze(-1)

X = torch.cat((X, X_batch), dim=0)
Y = torch.cat((Y, Y_batch), dim=0)

print(f'Evaluations: {len(Y)}, Minimum: {-Y.max()} at {unnormalize(X[Y.argmax()], bounds)}')
```

```python
plot_predictions_cut(X, Y, model, at=0.5, size=0.1)
```

See how for just 5 dimensions, even with a pretty large cut, we barely have any points to show? Every additional dimension makes the optimization problem drastically more difficult. 

Finally, we can plot the convergence of the algorithm. Let us define a new function for this.

The maximum value is `0` at `x = [0, 0, ..., 0]`. You should run at least 200 iterations to observe the convergence. See how it struggles to improve after a few iterations?

```python
def plot_convergence(Y):
    plt.figure()
    plt.scatter(torch.arange(len(Y)), Y, label="Samples")
    plt.scatter(Y.argmax(), Y.max(), label="Maximum")
    plt.xlabel('Evaluations')
    plt.ylabel('Value')
    plt.yscale("symlog")
    plt.legend()
```

```python
plot_convergence(Y)
```

# Turbo

See how for larger dimensions (dim >= 5), the algorithm really struggles to improve after a few interations? In such a large space, there will always be some areas where uncertainty is high and ExpectedImprovement will tend to do more exploration in these areas rather than exploiting the predicted mean. In other words, when working with many dimensions, ExpectedImprovement does more exploration than exploitation and fails to improve the bets value.

To solve this problem, we can use the Turbo algorithm (see: https://botorch.org/tutorials/turbo_1)

**Turbo is simple: every time we fail to improve the best value 3 iterations in a row, we *halve* the search space.**

**Note: this guarantees convergence, however there is a higher possibility that we get stuck in a local minima. But we can tune the convergance by changing the `failure_tolerance` parameter, which is `3` by default.**


## Turbo search space

The restricted search window is further weighted by the lengthscales of the kernel. We define a function that returns the new search bounds given the length of the Turbo window and the lengthscales of the kernel.

```python
def get_turbo_bounds(X, Y, model, length):
    # Scale the TR to be proportional to the lengthscales
    x_center = X[Y.argmax(), :].clone()

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
```

## Example

Here once again we look at the Ackley function in multiple dimensions.

```python
# Init samples
dim = 5
n_init = 100
batch_size = 5

# Turbo parameters
failure_counter = 0
failure_tolerance = 3 # setting this to torch.inf will never shrink the search range
turbo_length = 1.0 # start with 1.0 to cover the whole unit-hypercube and progressively shrink it

# Best value so far
Y_max: float = -torch.inf

sobol = SobolEngine(dimension=dim, scramble=True)

bounds = [-4, 4]
ackley_function = Ackley(dim=dim, negate=True)

X = sobol.draw(n=n_init) # the data is alreay in the correct shape, no need to unsqueeze
Y = ackley_function(unnormalize(X, bounds)).unsqueeze(-1)
```

```python
# Loop (ctrl+enter to run the same cell)

# Increase failure_counter if no improvement to the lossfunction
if Y.max().item() > Y_max: 
    failure_counter = 0
    Y_max = Y.max().item()
else:
    failure_counter += 1

# Shrink Turbo search range when failure_counter = failure_tolerance
if failure_counter == failure_tolerance:
    turbo_length /= 2.0
    failure_counter = 0

model = fit_SingleTaskGP(X, Y)

# Calculate Turbo search range using the covariance matrix
tr_lb, tr_ub = get_turbo_bounds(X, Y, model, turbo_length)

X_batch = generate_batch_ExpectedImprovement(X, Y, model, batch_size, tr_lb, tr_ub)
Y_batch = ackley_function(unnormalize(X_batch, bounds)).unsqueeze(-1)

X = torch.cat((X, X_batch), dim=0)
Y = torch.cat((Y, Y_batch), dim=0)

print(f'Evaluations: {len(Y)}, Failures: {failure_counter}, Turbo length: {turbo_length}, Minimum: {-Y.max()}')
```

```python
plot_predictions_cut(X, Y, model, at=0.5, size=0.01)
```

Finally, we can plot the convergence of the algorithm. The maximum value is `0` at `x = [0, 0, ..., 0]`.

```python
plot_convergence(Y)
```

# Turbo + Thompson sampling

Thompson Sampling is another method for generating new batches of points. It is particularly useful with many dimensions because it can generate many points quickly. It works well even when the search window gets very small and when the surrogate model is not so great, which happens ometimes when most of the samples are concentrated in a small areas. Thomspon sampling does more exploitation than exploration by design, which is what we want here with many dimensions. This means we also have a higher probability of getting stuck in a local minima. Therefore when using Thompson sampling it is recommended to increase the number of init samples.


## Thompson sampling

Let us compare the convergence for the same number of iterations.

```python
from botorch.generation import MaxPosteriorSampling

def generate_batch_ThompsonSampling(X, Y, model, batch_size, tr_lb, tr_ub):
    n_candidates = 2000

    dim = X.shape[-1]

    sobol = SobolEngine(dim, scramble=True)
    pert = sobol.draw(n_candidates)
    pert = tr_lb + (tr_ub - tr_lb) * pert

    # Create a perturbation mask
    prob_perturb = min(20.0 / dim, 1.0)
    mask = torch.rand(n_candidates, dim) <= prob_perturb
    ind = torch.where(mask.sum(dim=1) == 0)[0]

    if len(ind) > 0:
        mask[ind, torch.randint(0, dim - 1, size=(len(ind),))] = 1
        
    # Create candidate points from the perturbations and the mask
    x_center = X[Y.argmax(), :].clone()
    X_cand = x_center.expand(n_candidates, dim).clone()
    X_cand[mask] = pert[mask]

    # Sample on the candidate points
    thompson_sampling = MaxPosteriorSampling(model=model, replacement=False)
    
    with torch.no_grad():  # We don't need gradients when using TS
        X_batch = thompson_sampling(X_cand, num_samples=batch_size)

    return X_batch
```

## Benchmark comparison

```python
def optimize(X, Y, f, bounds, n_init, batch_size, batch_generator, failure_counter, failure_tolerance, turbo_length, Y_max, n):
    while len(X) < n:
        # Increase failure_counter if no improvement to the lossfunction
        if Y.max().item() > Y_max: 
            failure_counter = 0
            Y_max = Y.max().item()
        else:
            failure_counter += 1
        
        # Shrink Turbo search range when failure_counter = failure_tolerance
        if failure_counter == failure_tolerance:
            turbo_length /= 2.0
            failure_counter = 0
        
        model = fit_SingleTaskGP(X, Y)
        
        # Calculate Turbo search range using the covariance matrix
        tr_lb, tr_ub = get_turbo_bounds(X, Y, model, turbo_length)
        
        X_batch = batch_generator(X, Y, model, batch_size, tr_lb, tr_ub)
        
        Y_batch = f(unnormalize(X_batch, bounds)).unsqueeze(-1)
        
        X = torch.cat((X, X_batch), dim=0)
        Y = torch.cat((Y, Y_batch), dim=0)
        
        print(f'Evaluations: {len(Y)}, Failures: {failure_counter}, Turbo length: {turbo_length}, Minimum: {-Y.max()}')# at {X[Y.argmax()]}')

        if turbo_length < 1e-6:
            break
            
    return X, Y
```

```python
def benchmark(f, dim, bounds, batch_generator, n):
    # Init samples
    n_init = 300
    batch_size = 10
    
    # Turbo parameters
    failure_counter = 0
    failure_tolerance = 3 # setting this to torch.inf will never shrink the search range
    turbo_length = 1.0 # start with 1.0 to cover the whole unit-hypercube and progressively shrink it
    
    # Best value so far
    Y_max: float = -torch.inf
    
    sobol = SobolEngine(dimension=dim, scramble=True)
    
    X = sobol.draw(n=n_init) # the data is alreay in the correct shape, no need to unsqueeze
    Y = f(unnormalize(X, bounds)).unsqueeze(-1)
    
    X, Y = optimize(X, Y, f, bounds, n_init, batch_size, batch_generator, failure_counter, failure_tolerance, turbo_length, Y_max, n)

    model = fit_SingleTaskGP(X, Y)

    return X, Y, model
```

```python
from botorch.test_functions.synthetic import Rosenbrock

dim = 20

# bounds = [-4, 4]
# f = Ackley(dim=dim, negate=True)

bounds = [-5,10]
f = Rosenbrock(dim=dim, negate=True)

X_ei, Y_ei, model_ei = benchmark(f, dim, bounds, generate_batch_ExpectedImprovement, 1000)
X_ts, Y_ts, model_ts = benchmark(f, dim, bounds, generate_batch_ThompsonSampling, 1000)
```

```python
plot_predictions_cut(X_ei, Y_ei, model_ei, at=0.5, size=0.1)
plot_predictions_cut(X_ts, Y_ts, model_ts, at=0.5, size=0.1)

plot_convergence(Y_ei)
plot_convergence(Y_ts)
```

```python

```
