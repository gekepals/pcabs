# pcabs

# Description of pcabs
This code adds information about a constructive tau-abstraction to the Markov Equivalence Class. 
This allows the code to direct edges that were undirected in the Markov Equivalence Class. 
The output of the code is the Markov Abstraction Equivalence Class.
The code removes all models on which a constructive tau-abstraction is not possible. This means that on all models that are left in the Markov Abstraction Equivalence Class, a constructive tau-abstraction is possible. 
See Geke Pals, *Markov Abstraction Equivalence Classes: Combining the Markov Equivalence Class and Constructive tau-Abstraction to restrict the number of Causal Models* (2019) for more details, provided in this folder as pdf. 

# How to use pcabs
With the following code, the Markov Abstraction Equivalence Class will be created.

You need to be sure that you have added the pcalg and pcabs packages to your project:
```
library(pcalg)
library(pcabs)
```

First, you add the original causal model ML. This is a causal model consisting of 5 variables: X1, X2, X3, X4, X5. The data of ML is simulated with the following code:

```
x1 <- rnorm(4000)
x2 <- x1 + rnorm(4000)
x3 <- x2 + rnorm(4000)
x4 <- x2 + rnorm(4000)
x5 <- x3 + x4 + rnorm(4000)
dat <- cbind(x1, x2, x3, x4, x5)
```

As input for pcabs, we need the pdag of ML. This is a tabular representation of the Markov Equivalence Class of ML, that stores which variables are connected with each other. We create the pdag of ML by first creating the skeleton of ML with the pcalg package. The pdag is then created with the find_pattern function from pcabs.

```
labels <- colnames(dat)
n <- nrow(dat)

skel <- skeleton(suffStat = list(C = cor(dat), n=n),
                 indepTest = gaussCItest, alpha = 0.01,
                 labels = labels, verbose = TRUE)

pdag <- find_pattern(skel)
```

We can plot this pdag, to see the Markov Equivalence Class of ML:

```
mec <- as(pdag, "graphNEL")
plot(mec)
```

Next, we define the abstraction-groups of M in a list. 
These are the clusters in M that are abstracted to high-level variables. 
These clusters are defined in the partition P of M.

```
a <- list("x1", "x2")
b <- list("x3", "x4", "x5")
abs_groups <- list(a, b)
```

Now, we can start pcabs. The final output of pcabs is new_pdag, which is the Markov Abstraction Equivalence Class of ML. We can plot this class as well.

```
new_pdag <- lead_abs(pdag, abs_groups)
pdag_plot <- as(new_pdag, "graphNEL")
plot(pdag_plot)
```
