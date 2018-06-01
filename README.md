# Rotterdam_code

Meg's code for testing reference free methods on GenR data, December 2017

| Code                         | Function                                                                | Input                                                                                            | Output                                              |   |
|------------------------------|-------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------|-----------------------------------------------------|---|
| Normalization_forGenR        | Normalization (bmiq)                                                    | Raw betas                                                                                        | BMIQ normalized betas                               |   |
| Initial_deconvolution_edited | Deconvolution of betas                                                  | Sorted Library, normalized raw betas                                                             | Louie_predicted_WB_celltypes_Oct26_fixed_code.rdata |   |
| FACS_decon_correct.R         | Correct the betas using the facs and deconvolution counts               | WB_betas_BMIQ_comabt_alloutliersremoved.rdata; GenR_cord_deconvolution_predicted_celltypes.rdata | method_corrected_betas.Rdata                        |   |
| Reference_free_methods.R     | Calculate the components and corrected based on ref free methods        | WB_betas_BMIQ_comabt_alloutliersremoved.rdata                                                    | counts.Rdata; adj.residuals_refactor.Rdata          |   |
| R2.R                         | Regress facs counts with all count types and look at explained variance | all count types                                                                                  | R-squared values from regressions                   |   |
| Error.R                      | MAE of betas against gold-standard FACS PCA corrected betas             | Corrected Betas                                                                                  | MAE values                                          |   |
| madcalc.R                    | Mad calcs between the count data                                        | all count types                                                                                  | mad values                                          |   |
| *_ewas.R                     | EWAS on either GA or sex                                                | All corrected betas                                                                              | p values for all CpGs                               |   |
