# bb-drones-pools-201911

![](graph.svg)

- Use honeybee-genotype-pipeline to genotype the pools and drones together
- pull the drones out of the filtered VCF
- create a "consensus read" for each drone
- use the consensus read with `whatshap phase` to phase the pool reads
