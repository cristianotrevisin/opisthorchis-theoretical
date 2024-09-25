# A spatially explicit model of the dynamics of *Opisthorchis viverrini* spread by C. Trevisin et al. (2024)

This repository contains the codes used to perform the analyses present in the aforementioned paper. 

## Scripts to run to generate the figures
- `experiment_1.m` carries out the **Preliminary comparison** (Figure 3 of the main text, Figures S1 and S2 of the Supporting Information. The reader will have to select on which OCN the analysis should be carried out (via the `ocnmap` toggle).
- `experiment_2.m` investigates the **Effect of fish trading** (Figure 4 of the main text, Figures S4 and S5 of the Supporting Information. The reader will have to select on which OCN the analysis should be carried out (via the `ocnmap` toggle).
- `experiment_3.m` investigates the **Effect of fish mobility via river corridors** (Figure 5 and 6 of the main text, Figures S6-S9 of the Supporting Information. The reader will have to select whether to carry out the analyses with the random allocation or downstream accumulation of the human population (via the `downstream_accumulation` boolean).
- `experiment_4.m` carries out the **Sensitivity analysis** (Figure 7 of the main text)
- `experiment_5.m` investigates the **Joint effect of fish market and mobility** (Figure 8 of the main text, Figures S10 and S11 of the Supporting Information. The reader will have to select on which OCN the analysis should be carried out (via the `ocnmap` toggle).

## Helper scripts
- `build_OCN.m` builds the geographical setup based on the OCN generated via the [OCNet Package](https://github.com/lucarraro/OCNet).
- `build_setup.m` builds the epidemiological setup based on whether the human population is randomly allocated to the nodes or downstream accumulated.
- `draw_OCN.m` makes the beautiful maps featured in the paper
- `drawborders.m` helps draw the borders among communities
- `assignSC.m` helps assign attributes to pixels so they can appropriately be coloured
- `compute_EE.m` finds the endemic equilibria for a system of disconnected nodes
- `calculateWRecursive.m` calculates the hydrological connectivity matrix (which guarantees a demographic equilibrium in the fish population)
- `sensitivity_analysis_EE.m`contains the instructions for the sensitivity analysis carried out via `experiment_4.m`
- `zipf.m` generates a Zipf-distributed human population
- `common_parameters.m` calls a dictionary containing the common parameters of the model.

## Additional remarks
- The reader may develop additional OCN maps by using the R script `generate_OCN.R` under `dataOCN`. 
