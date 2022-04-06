scVI_arg_list <- list(n_layers = c(2,3,4),
                      n_hidden = c(128,256),
                      n_latent = c(50,65,80,95,110),
                      dropout_rate =c(0.01,0.05,0.1),
                      max_epochs=c(20,50,100,150,200,300,400),
                      dispersion = "gene",
                      gene_likelihood = "zinb",
                      early_stopping = c(FALSE))

# save
scUtils::writeList_to_JSON(scVI_arg_list,filename = "data/parameters_scvi_1.json")
