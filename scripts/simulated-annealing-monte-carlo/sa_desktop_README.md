### Step 3: Run the Simulated Annealing Monte Carlo Sampling

After completing the previous steps, you should have generated the metadata file:

```text
generated_barcodes.tsv
```

Alternatively, you can use the pre-generated file included in the repository.

This metadata file is the only input required to run the simulated annealing algorithm.

Run:

```text
scripts/simulated-annealing-monte-carlo/sa_desktop.py
```

Before running the script, update the input file path:

```python
tsv_file = "/path/to/generated_barcodes.tsv"
```

The main simulation parameters can be modified in the following section of the script:

```python
if __name__ == "__main__":
    run_spatial_simulated_annealing(
        tsv_file="/path/to/generated_barcodes.tsv",
        num_runs=1068,
        k_neighbors=150,
        spatial_weight=0.6,
        energy_weight=0.4,
        outer_iterations=200,
        inner_iterations=50,
    )
```

### Adjustable Parameters

#### Number of simulation runs

```python
num_runs = 1068
```

This parameter specifies how many simulated annealing trajectories are performed.

In our study, `1068` corresponds to the number of **OPC (Oligodendrocyte Precursor Cell)** cells. Because OPCs are multipotent progenitor cells that can differentiate into **ASC** and **OGC**, the simulated annealing algorithm is initialized once from each OPC cell. Consequently, one optimization trajectory is generated for every OPC cell.

#### Number of nearest neighbors

```python
k_neighbors = 150
```

This parameter determines the number of nearest neighbors considered in the MDS embedding when evaluating candidate transitions.

Changing this value alters the neighborhood structure and therefore changes the simulated annealing trajectories and the resulting figures.

#### Spatial and energy weights

```python
spatial_weight = 0.6
energy_weight = 0.4
```

These parameters control the relative importance of the spatial proximity and energy terms in the objective function.

Increasing a weight increases the influence of the corresponding term during optimization.

#### Annealing schedule

```python
outer_iterations = 200
inner_iterations = 50
```

* `outer_iterations` specifies the number of temperature levels in the annealing schedule (i.e., how many times the temperature is decreased).
* `inner_iterations` specifies the number of optimization steps performed at each temperature.

Increasing either parameter generally results in a more thorough exploration of the solution space but also increases the computational time.

#### Temperature parameters

The simulated annealing schedule is also controlled by the following parameters:

```python
T_start = 0.2
T_min = 1e-6
alpha = 0.98
```

* `T_start` is the initial temperature.
* `T_min` is the minimum temperature at which the algorithm terminates.
* `alpha` is the cooling factor applied after each outer iteration.



> **Note:** Changing any of these parameters modifies the behavior of the simulated annealing algorithm. As a result, the generated trajectories, energy evolution, and reproduced figures may differ from those reported in the paper, although the overall methodology remains the same.

