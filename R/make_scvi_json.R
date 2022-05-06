scVI_arg_list <- list(n_layers = c(3,4),
                      n_hidden = c(256), # 128
                      n_latent = c(65),# ---> also defined in main parameter json! but i am overwriting atm
                      dropout_rate =c(0.1),
                      max_epochs=c(300,400),
                      dispersion = "gene",
                      gene_likelihood = "zinb",
                      early_stopping = c(FALSE))

# save
scUtils::writeList_to_JSON(scVI_arg_list,filename = "data/parameters_scvi_neurons_2.json")
