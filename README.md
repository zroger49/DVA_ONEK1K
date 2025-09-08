# DVA ONEK1K


Gene expression variability across cells is increasingly recognized as a fundamental dimension of transcriptomic regulation, capturing cell-to-cell heterogeneity beyond changes in average expression levels. Despite advances in single-cell RNA sequencing (scRNA-seq), quantifying and comparing expression variability remains challenging due to data sparsity, zero inflation, and variation in sequencing depth and capture efficiency. Additionally, strong mean–variance dependency can obscure true variability effects and generate spurious associations driven by changes in mean expression.

To address these issues, the **Mapping genetic Effects On inTra-Individual Variability of gene Expression (MEOTIVE)** framework estimates gene-level dispersion using a Cox–Reid adjusted maximum likelihood approach, providing a robust proxy for intra-individual variability. Importantly, this measure is independent of mean expression, yielding an unbiased estimator. MEOTIVE has been used to uncover genetic variants associated with expression variability, highlighting their role in regulatory plasticity—particularly in immune cells where stochastic transcriptional dynamics are functionally relevant.

---

## Project Overview

In this work, the MEOTIVE framework was adapted to investigate **differential expression variability** in the context of **aging** and **sex differences**.  

### Workflow
1. **Dispersion Estimation**  
   - Gene-level dispersion was estimated per donor for each cell population.
   - Generalized linear models (GLMs) was used to model the dispersion across individuals.
   - A mean expression filter was applied to remove lowly expressed genes, and p-values were readjusted. The threshold was derived from simulations accounting for the number of cells per donor and the number of donors.  

2. **Group Comparisons**  
   - Comparison between male vs female, age, and age stratified by male and female.
   - Confounders were adjusted for in the models.
     
3. **Robustness**  
   - Results were validated with both individual-level and cell-level permutations.  

## Limitations

- Dispersion estimates are inflated for lowly expressed genes.  
- Filtering thresholds increase for smaller cell populations.  
- Current dataset limited the analysis to highly expressed genes in abundant cell types.  


---

## License

This project is made available under the MIT License.
