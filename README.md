# Neural field models for latent state inference: Application to large-scale neuronal recordings

Rule, M. E., Schnoerr, D., Hennig, M. H., & Sanguinetti, G. (2019). Neural field models for latent state inference: Application to large-scale neuronal recordings. [PLOS Computational Biology](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007442), 15(11), e1007442. [doi:  https://doi.org/10.1371/journal.pcbi.1007442]( https://doi.org/10.1371/journal.pcbi.1007442).

### Abstract

Large-scale neural recording methods now allow us to observe large populations of identified single neurons simultaneously, opening a window into neural population dynamics in living organisms. However, distilling such large-scale recordings to build theories of emergent collective dynamics remains a fundamental statistical challenge. The neural field models of Wilson, Cowan, and colleagues remain the mainstay of mathematical population modeling owing to their interpretable, mechanistic parameters and amenability to mathematical analysis. Inspired by recent advances in biochemical modeling, we develop a method based on moment closure to interpret neural field models as latent state-space point-process models, making them amenable to statistical inference. With this approach we can infer the intrinsic states of neurons, such as active and refractory, solely from spiking activity in large populations. After validating this approach with synthetic data, we apply it to high-density recordings of spiking activity in the developing mouse retina. This confirms the essential role of a long lasting refractory state in shaping spatiotemporal properties of neonatal retinal waves. This conceptual and methodological advance opens up new theoretical connections between mathematical theory and point-process state-space models in neural data analysis.

### Summary

Developing statistical tools to connect single-neuron activity to emergent collective dynamics is vital for building interpretable models of neural activity. Neural field models relate single-neuron activity to emergent collective dynamics in neural populations, but integrating them with data remains challenging. Recently, latent state-space models have emerged as a powerful tool for constructing phenomenological models of neural population activity. The advent of high-density multi-electrode array recordings now enables us to examine large-scale collective neural activity. We show that classical neural field approaches can yield latent state-space equations and demonstrate that this enables inference of the intrinsic states of neurons from recorded spike trains in large populations.

### Organization

This repository contains a Matlab implementation and demonstration of the methods described in *Neural field models for latent state inference: Application to large-scale neuronal recordings.*

 - `./NFCP` Matlab source code
 - `./Examples` Demonstration scripts
 - `./Documentation` HTML files documenting the Matlab files
 - `./Gittools` utility scripts for github
 
