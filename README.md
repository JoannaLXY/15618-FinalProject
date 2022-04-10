# Barnes Hut algorithm for N-Body problem
## Summary
We are going to implement Barnes Hut algorithm for N-Body problem using OpenMP, and perform a detailed analysis on sequential and parallel versions.


## Background
In physics, the N-body problem is to predict the individual movements of a group of stars and planets that are influenced by each other through gravitational forces. The gravitational forces change the speeds and positions of the objects, and in turn, the positions will change the forces among them. Our goal is that given the initial information about these objects (position, velocity and time), predict its true orbits and trajectory for any future time. 

Barnes Hut algorithm consists of two parts that can be parallelized: tree building and force computing. In tree building, it recursively divides the n-bodies into groups and build the tree. In the force computing, it traverses the tree to calculate the forces and update the positions. Both parts could be benefit from the parallelism since the computation for each body is independent of others. We use the pseudocode provided in lecture as the basic idea.

```
for each time step in simulation:
    build tree structure
    compute (aggregate mass, center-of-mass) for interior nodes 
    for each particle:
        traverse tree to accumulate gravitational forces update particle position based on gravitational forces
```

## Challenge

The naive approach to this problem is very straightforward: for each object, calculate forces that other objects have on it, and use the combined force to calculate the future orbits. However, doing so results in a O(N^2) time complexity as all particles interact with all others. Hence, we want to use Barnes Hut algorithm to improve the performance without comprising the precision too much. However, Barnes Hut algorithm also imposes many difficulties as follows:
1. The workload distribution is not uniform. Some part with high density may need more computation and communication.
2. The workload for each body is not the same. So workload assignment based on number of bodies may result in unbalanced workload for each processor.
3. The workload distribution is not static over time. The bodies would move so the costs and communication patterns would change over time.

## Resources
We will use GHC and PSC machines to develop our project. We will start from scratch, but use this website[1] as a reference for logic. We also use the test files from [this GitHub repository](https://github.com/chindesaurus/BarnesHut-N-Body/tree/master/inputs) for testing and benchmark.

## Goals and Deliverables

PLAN TO ACHIEVE ("100\%"): We plan to finish implementing the sequential version and OpenMP version of the algorithm. Also, we will compare the performance between sequential and parallel version.

Extra goal or two ("125\%"): Besides the sequential and OpenMP version, we could implement the CUDA version if we are ahead of the schedule. 

Goals in case the work goes more slowly ("75\%"): We finish the sequential version and implement a workable parallel solution which may not speed up as much as expected. 

As for demo, we plan to visualize our outputs which is the evolution of particle positions in video. We will also show the speedup graphs compared with sequential version.

If we have time, we can construct the tree parallel using Morton Key method (https://github.com/scbrubaker02/hpc-barnes-hut)

## Platform Choice
We use GHC for developing and simple testing. We choose GHC because GeForce RTX 2080 comes with the machines, on top of which we can build our CUDA version. We also choose PSC because it can provide many nodes, which helps up exploit the power of parallelism.

## Schedule
Week of 3.28 Finish sequential version 

Week of 4.4 Produce workable OpenMP version and milestone report

Week of 4.11 Improve OpenMP version

Week of 4.18 Visualization

Week of 4.25 Analysis, benchmarking and final report

## MILESTONE

### Progress

We have finished the sequential version of Barnes-Hut algorithm and are still working on the parallel implementation using OpenMP. We implemented algorithm using the tree-based method according to this website (https://beltoforion.de/en/barnes-hut-galaxy-simulator/) in C++.

Besides, when calculating the particle movements, we also save the status of each particle in csv file, and then use visualization code to produce an animation.


### How we are doing with respect to the goals
We plan to implement three versions of barnes-hut algorithm: sequential, parallel with OpenMP and parallel with CUDA. We have finished the implementation of the sequential version and are still working on the openMP version. We believe we can produce all the deliverables we have written in our proposal.


### The new list of the goals is:
Continue on openMP implementation

Finish CUDA implementation

Finish poster

Finish final report


### Poster session
We plan to provide performance analysis for different parallel methods and parameters using graphs and figures during the poster session. We would also provide a live demo to produce animations of some sample n-body movements. 


### Concerns
The concerns we have so far are:
1. Correctness: Currently we don't have a generic method of evaluating implementation correctness. For some very simple tests, we can calculate the result by hand, but this is very inefficient and unrealistic. We can have a rough idea from the animation, but still can not determine whether the implementation is correct or not.
2. Parallelism: Currently we plan to distribute the force computation to different working threads. However, as the number of particles increases, the tree building phase will take on a larger portion, and we need to figure out how to speed up this part.

### Updated Schedule
Week of 4.11 Finish and improve OpenMP version

Week of 4.18 Finish and improve CUDA version

Week of 4.25 Analysis, benchmarking and final poster & report


## Reference
[1] The barnes-hut galaxy simulator. https://beltoforion.de/en/barnes-hut-galaxy-simulator/. Accessed: 2020-03-21.
