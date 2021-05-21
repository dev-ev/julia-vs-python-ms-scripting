# julia-vs-python-ms-scripting
## Performance of Julia and Python for reading and matching mass spectrometry data


The example MGF data file is provided in this repository, it contains 1000 tandem mass spectra from a Baker's yeast sample. The MS data is taken from the [project PXD021218 on PRIDE repository website](https://www.ebi.ac.uk/pride/archive/projects/PXD021218/).

Reading the MGF file with the Python function:
     
    fname = 'Yeast_1000spectra.mgf'
    %%timeit -r 20
    load_mgf(fname)
  
    290 ms ± 7.45 ms per loop (mean ± std. dev. of 20 runs, 1 loop each)
  
Implementing the same logic in Julia:

    @benchmark loadmgf(fname)
  
    BenchmarkTools.Trial: 
      memory estimate:  48.80 MiB
      allocs estimate:  911621
      --------------
      minimum time:     108.247 ms (6.60% GC)
      median time:      112.938 ms (5.28% GC)
      mean time:        114.375 ms (5.78% GC)
      maximum time:     126.844 ms (4.72% GC)
      --------------
      samples:          44
      evals/sample:     1
  
Matching the spectral data against the set of mass deltas, implemented with Python/numpy:

     %%timeit -r 5
     find_matches(res, singleResDeltas, 1e-5)
  
     2.49 s ± 135 ms per loop (mean ± std. dev. of 5 runs, 1 loop each)
  
Implementing the same matching logic as closely as possible in Julia:

     @benchmark matches(res, single_res_Δ, 1e-5) samples=5
  
     BenchmarkTools.Trial: 
         memory estimate:  5.84 GiB
         allocs estimate:  34240
         --------------
         minimum time:     1.604 s (8.30% GC)
         median time:      1.787 s (7.45% GC)
         mean time:        1.734 s (7.71% GC)
         maximum time:     1.812 s (7.14% GC)
         --------------
         samples:          3
         evals/sample:     1
    
Julia allows us to fuse several matrix operations into one. Applying the modified matching algorithm with the fused matrix operation in Julia:
  
     @benchmark matchesfuse(res, single_res_Δ, 1e-5) samples=5
  
     BenchmarkTools.Trial: 
         memory estimate:  1.54 GiB
         allocs estimate:  34246
         --------------
         minimum time:     1.036 s (3.13% GC)
         median time:      1.078 s (3.21% GC)
         mean time:        1.071 s (3.31% GC)
         maximum time:     1.087 s (3.72% GC)
         --------------
         samples:          5
         evals/sample:     1
