# Sparse-PCA-MLSE
Sparse PCA on white matter connectome data for MLSE project

For my postdoc project using the MLSE data, I used connectomics to find patterns of white matter connectivity degeneration in each patient group. Given the number of connections, I used sparse PCA from the PMA package in R to find meaningful patterns of change in these groups.

Sparse PCA requires hyper parameter (sparsity parameter) tuning. I found this resource - Tibshirani group (the one that initially proposed sparse PCA) have a package of their own called PMA. It is described in more detail here: https://academic.oup.com/biostatistics/article/10/3/515/293026. That package has a sparse PCA and a cross-validation algorithm that helps with tuning parameter selection. The package also allows to do orthogonal PCA if needed. As part of the PMA package, the SPC algorithm runs a sparse PCA with penalty applied to columns (not rows)

The excerpt of the code for sparse PCA is here:
Main steps involve: (i) choose number of components = 36 as we derived from component selection (code not shown here); (ii) for sparsity parameter tuning, run SPC.CV; (iii) store the sumabsv sparsity parameter as a variable; (iv) plug it into the SPC; (v) save the u output as rotated SCORES, the v output as rotated LOADINGS, and the prop.var.explained as VAR EXPL

load the library<br />
`library(PMA)`

First, I select a sparsity parameter. This is represented by the sum of absolute values (sumabsv). The default is the sqrt(ncol(x)), in this case, square root of ncol of our df, which is 81. If I choose this value, then the solution is the least sparse one. The lower the sumabsv value, the more sparse the solution is. To find the best sumabsv, I tell the cross-validation algorithm to consider everything between 1 and the sqrt(ncol(matrix)) value. The idea is that the K has already been selected through our other venetian blind cv method.

`out.cv.x<-SPC.cv(as.matrix(df), # input the matrix`<br />
                `sumabsvs=seq(1, sqrt(ncol(df)), # for hyperparameter tuning, set the limits between 1 and sqrt of ncols (which is the highest you can have - see help)`<br />
              `niter=100, nfolds=10, orth = T) # 100 iterations, 10 folds, and  components to be orthogonal`<br />

`out.cv.x$bestsumabsv`<br />

The suggested sumabsv value is 73. This is close enough to the sqrt(ncol(mlse_conn_df_FA_85_resid_forPCA_withSite.resid)), which is 81, so lets plug the `out.cv.x$bestsumabsv` as the tuning parameter in to the SPC solution.<br />
`out.tmp<-SPC(as.matrix(df), sumabsv = out.cv.x$bestsumabsv, K=36, orth = T,niter=100, center = F)`<br />

`out.tmp$prop.var.explained`<br />
91% variance explained

Plot the variance explained. Note, that I have done sparse PCA for 5 different indices (FA, SC, AD, RD, ADC) from the connectome, concatenated and reshaped the data into a new df called df_melt (code not shown here)<br />

`df_sparse_varexp <- ggplot(df_melt, aes(x=as.factor(Component), y=value, group=VarExpIndex)) +` <br />
  `geom_step(aes(colour=VarExpIndex)) +`<br />
  `geom_point(aes(colour=VarExpIndex)) +`<br />
  `geom_label_repel(force_pull   = 0, # do not pull toward data points`<br />
                    `nudge_x      = 0.05,`<br />
                    `direction    = "x",`<br />
                    `angle        = 90,`<br />
                  ` hjust        = 1,`<br />
                   `segment.size = 1,`<br />
                   `min.segment.length = 0,`<br />
                   `lineheight=5,`<br />
                   `xlim=c(NA, NA),`<br />
                  `aes(label=round(value, digits=0)), `<br />
             `data = . %>%`<br />
              `filter(row_number() %% 1 == 0 &`<br />
                 `row_number() %% 4 == 0)) +`<br /> 
  `geom_label_repel(force_pull   = 0,`<br /> 
                   `nudge_x      = 0.05,`<br />
                   `direction    = "y",`<br />
                   `angle        = 90,`<br />
                   `hjust        = 0,`<br />
                   `segment.size = 1,`<br />
                   `min.segment.length = 0, `<br />
                   `data=subset(df_melt1, Component==1 & VarExpIndex=='Component variance'),`<br />
                   `aes(label=round(value, digits=0))) +`<br />
   `scale_y_continuous(breaks=seq(0, 100, 10), limits=c(0, 100)) +`<br />
  `facet_wrap(~ConnIndex, scales="free") +`<br />
  `xlab("Component") +`<br />
  `ylab("Variance Explained") +`<br />
  `theme_bw() +`<br />
  `theme(strip.text.x = element_text(size=11),`<br />
        `legend.background = element_blank(),`<br />
        `legend.box.background = element_rect(colour = "black")) +`<br />
  `coord_flip() +`<br />
  `labs(color="Variance\nExplained")`<br />
  
  ![SparsePCA_VarExplain_85Percent](https://user-images.githubusercontent.com/88196987/232824072-6da04e30-df68-4268-8038-a1f841a20054.jpeg)

Show the plot for sparse PCA loadings for all 5 indices

![LoadingDistribution](https://user-images.githubusercontent.com/88196987/232826234-2f39d4e3-41bd-4b80-81d0-4bb0a4476fe4.jpeg)


Next, I conduct regressions between sparse PCA loadings (connectome) and regular PCA loadings (behaviour) to find associations between these two datasets (code not shown) and plot these results in a tile plot (plot below shows only one PCA factor).
![Comp3](https://user-images.githubusercontent.com/88196987/232825772-94c5f075-c0b3-448e-a394-b05c01a221d2.jpeg)

Project the same plot from above but into the brain (see right panel for white matter)

![Comp2](https://user-images.githubusercontent.com/88196987/232826504-655f6a44-11eb-4ac2-bfd9-caec095617a3.JPG)
