# Sparse-PCA-MLSE
Sparse PCA on white matter connectome data for MLSE project

For my postdoc project using the MLSE data, I used connectomics to find patterns of white matter connectivity degeneration in each patient group. Given the number of connections, I used sparse PCA from the PMA package in R to find meaningful patterns of change in these groups.

The excerpt of the code for sparse PCA is here:
# SParse PCA requires hyper parameter (sparsity parameter) tuning
# I found this resource - Tibshirani group (the one that initially proposed sparse PCA) have a package of their own called PMA...
# ... it is described in more detail here: https://academic.oup.com/biostatistics/article/10/3/515/293026
# ... that package has a sparse PCA and a cross-validation algorithm that helps with tuning parameter selection...
# ... the package also allows to do orthogonal PCA if needed.
# as part of the PMA package, the SPC algorithm runs a sparse PCA with penalty applied to columns (not rows)

# finalise Sparse PCA Pipeline = 
# choose number of components = 36 as we derived from component selection (code not shown here)
# ... for sparsity parameter tuning, run SPC.CV...
# ... store the sumabsv sparsity parameter as a variable...
# ... plug it into the SPC...
# ... save the u output as rotated SCORES...
# ... save the v output as rotated LOADINGS...
# ... save the prop.var.explained as VAR EXPL

# load the library
library(PMA)

# First, I select a sparsity parameter...
# this is represented by the sum of absolute values (sumabsv)..
#... the default is the sqrt(ncol(x)), in this case, square root of ncol of our df, which is 81..
# ... if I choose this value, then the solution is the least sparse one...
# ... the lower the sumabsv value, the more sparse the solution is...

# to find the best sumabsv, I tell the cross-validation algorithm to consider everything between 1 and the sqrt(ncol(matrix)) value
# ... the idea is that the K has already been selected through our other venetian blind cv method...

out.cv.x<-SPC.cv(as.matrix(df), # input the matrix
                 sumabsvs=seq(1, sqrt(ncol(df)), # for hyperparameter tuning, set the limits between 1 and sqrt of ncols (which is the highest you can have - see help)
              niter=100, nfolds=10, orth = T) # 100 iterations, 10 folds, and  components to be orthogonal

out.cv.x$bestsumabsv

# the suggested sumabsv value is 73...
# ... this is close enough to the sqrt(ncol(mlse_conn_df_FA_85_resid_forPCA_withSite.resid)), which is 81
# so lets plug the out.cv.x$bestsumabsv as the tuning parameter in to the SPC solution
out.tmp<-SPC(as.matrix(df), 
         sumabsv = out.cv.x$bestsumabsv,
         K=36, orth = T,
         niter=100, center = F)

out.tmp$prop.var.explained
# 91% variance explained

# plot the variance explained
# note, that I have done sparse PCA for 5 different indices (FA, SC, AD, RD, ADC) from the connectome, concatenated and reshaped the data into a new df called df_melt (code not shown here)

df_sparse_varexp<-ggplot(df_melt, aes(x=as.factor(Component), y=value, group=VarExpIndex)) +
  geom_step(aes(colour=VarExpIndex)) +
  geom_point(aes(colour=VarExpIndex)) +
   geom_label_repel(force_pull   = 0, # do not pull toward data points
                    nudge_x      = 0.05,
                    direction    = "x",
                    angle        = 90,
                   hjust        = 1,
                   segment.size = 1,
                   min.segment.length = 0,
                   lineheight=5,
                   xlim=c(NA, NA),
                   # max.iter = 1e4, max.time = 1,
   # min.segment.length = 0, box.padding = 0.5,
                  aes(label=round(value, digits=0)), # round the var exp val
             data = . %>%
               filter(row_number() %% 1 == 0 &
                 row_number() %% 4 == 0)) + # label every 4th row. if you want | row_number() == 0)) +
  geom_label_repel(force_pull   = 0, # repeat step again but this time just for the first component
                   nudge_x      = 0.05,
                   direction    = "y",
                   angle        = 90,
                   hjust        = 0,
                   segment.size = 1,
                   min.segment.length = 0, 
                   data=subset(df_melt1, Component==1 & VarExpIndex=='Component variance'), # ensures the number doesnt double up for both component & cumulative variance
                   aes(label=round(value, digits=0))) +
   scale_y_continuous(breaks=seq(0, 100, 10), limits=c(0, 100)) +
  facet_wrap(~ConnIndex, scales="free") +
  xlab("Component") +
  ylab("Variance Explained") +
  theme_bw() +
  theme(strip.text.x = element_text(size=11),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  #scale_x_discrete(limits = rev) +  #scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) # if you want to skip intermediate labels  #scale_x_discrete(guide=guide_axis(n.dodge = 2)) # if you want every alternate label to be one line below
  coord_flip() +
  labs(color="Variance\nExplained")
  
  ![SparsePCA_VarExplain_85Percent](https://user-images.githubusercontent.com/88196987/232824072-6da04e30-df68-4268-8038-a1f841a20054.jpeg)

# show the plot for sparse PCA loadings for all 5 indices

![LoadingDistribution](https://user-images.githubusercontent.com/88196987/232826234-2f39d4e3-41bd-4b80-81d0-4bb0a4476fe4.jpeg)


# Next, I conduct regressions between sparse PCA loadings (connectome) and regular PCA loadings (behaviour) to find associations between these two datasets (code not shown) and plot these results in a tile plot (plot below shows only one PCA factor).
![Comp3](https://user-images.githubusercontent.com/88196987/232825772-94c5f075-c0b3-448e-a394-b05c01a221d2.jpeg)

# project the same plot from above but into the brain (see right panel for white matter)

![Comp2](https://user-images.githubusercontent.com/88196987/232826504-655f6a44-11eb-4ac2-bfd9-caec095617a3.JPG)
