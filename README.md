# pcabs

# Description of pcabs
This code adds information about a constructive tau-abstraction to the Markov Equivalence Class. 
This allows the code to direct edges that were undirected in the Markov Equivalence Class. 
The output of the code is the pcabs-class: a class that lies between the Markov Equivalence Class and the Markov tau-Abstraction Equivalence Class. 
The code only removes models on which a constructive tau-abstraction is not possible. 
However, it cannot be guaranteed that a constructive tau-abstraction is possible on all models that are left in the pcabs-class.

# How to use pcabs
With the following code, the pcabs-Class will be created.

First, you add the original causal model M:

```
x1 <- rnorm(4000)
x2 <- x1 + rnorm(4000)
x3 <- x2 + rnorm(4000)
x4 <- x2 + rnorm(4000)
x5 <- x3 + x4 + rnorm(4000)
dat <- cbind(x1, x2, x3, x4, x5)
```

For this causal model, you can find the Markov Equivalence Class using the pcalg package. This is not necessary for pcabs.

```
labels <- colnames(dat)
n <- nrow(dat)

pc.fit <- pc(suffStat = list(C = cor(dat), n = n),
             indepTest = gaussCItest, alpha = 0.01,
             labels = labels, verbose = TRUE)
pc.fit
plot(pc.fit)
```

pcabs needs the skeleton of M to generate the pdag of M, which is a tabular representation of the edges and variables of M.

```
labels <- colnames(dat)
n <- nrow(dat)

skel <- skeleton(suffStat = list(C = cor(dat), n=n),
                 indepTest = gaussCItest, alpha = 0.01,
                 labels = labels, verbose = TRUE)
plot(skel)
```

With the skeleton, the pdag of M can be generated:

```
pdag <- find_pattern(skel)
```

Next, we define the abstraction-groups of M in a list. 
These are the clusters in M that are abstracted to high-level variables. 
These clusters are defined in the partition P of M.

```
a <- list("x1", "x2")
b <- list("x3", "x4", "x5")
abs_groups <- list(a, b)
```

Now, we can start pcabs. 

```
new_pdag <- lead_abs(pdag, abs_groups)
pdag_plot <- as(new_pdag, "graphNEL")
plot(pdag_plot)
```
