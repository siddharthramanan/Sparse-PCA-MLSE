# Sparse-PCA-MLSE

In this project, I use an in-house dataset containing measures of language performance and white matter brain connectivity changes in three patinet groups. As the data is high-dimensional, I use connectomics to find patterns of white matter connectivity degeneration in each patient group. Given the number of connections, I used sparse PCA from the PMA package in R to find meaningful patterns of change in these groups.<br />
This is an example picture of my workflow:
![MLSE_WM_Workflow](https://user-images.githubusercontent.com/88196987/232834166-b5b1d69a-2488-41ed-b4ab-ac28436a61b3.jpg)

Sparse PCA requires hyper parameter tuning. I use the PMA package in R from the the Tibshirani group which is described in more detail here: https://academic.oup.com/biostatistics/article/10/3/515/293026. PMA package has a sparse PCA and an accompanying cross-validation algorithm for parameter tuning. The package also allows to do orthogonal PCA if needed.

An excerpt of my sparse PCA code pipeline is here:
The main steps involves: (i) choose 36 components (this number was derived from an independent component selection algorithm not shown here); (ii) running SPC.cv command for sparsity parameter tuning; (iii) storing the sumabsv sparsity parameter as a variable; (iv) plug it into the sparse PCA; (v) save the u output as rotated SCORES, the v output as rotated LOADINGS, and the prop.var.explained as VAR EXPL

load the library<br />
`library(PMA)`

Selecting the sparsity parameter: In SPC.cv, this is represented by the sum of absolute values (sumabsv). The default is the sqrt(ncol(x)) which leads to the least sparse solution. The lower the sumabsv value, the more sparse the solution is. To find the best sumabsv, I tell the cross-validation algorithm to consider everything between 1 and the sqrt(ncol(matrix)) value.

`out.cv.x<-SPC.cv(as.matrix(df), # input the matrix`<br />
                `sumabsvs=seq(1, sqrt(ncol(df)), # for hyperparameter tuning, find the value in a sequence between limits between 1 and sqrt of ncols`<br />
              `niter=100, nfolds=10, orth = T) # 100 iterations, 10 folds, and  components to be orthogonal`<br />

`out.cv.x$bestsumabsv`<br />

The suggested sumabsv value is 73. I plug this as the tuning parameter in to the SPC solution.<br />
`out.tmp<-SPC(as.matrix(df), sumabsv = out.cv.x$bestsumabsv, K=36, orth = T,niter=100, center = F)`<br />

`out.tmp$prop.var.explained`<br />
91% variance explained

Plot the variance explained: Note, that I have done sparse PCA for 5 different indices (FA, SC, AD, RD, ADC) from the connectome, concatenated and reshaped the data into a new df called df_melt (code not shown here)<br />

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

I repeat the sparse PCA for all 5 indices and display the loadings in this plot below (plot code not displayed here).

![LoadingDistribution](https://user-images.githubusercontent.com/88196987/232826234-2f39d4e3-41bd-4b80-81d0-4bb0a4476fe4.jpeg)

Next, I choose the highest loading sparse PCs and use regressions to examine associations between these high loading sparse PCs and the regular PCA done on behavioural data (code not shown here). I plot these results in a tile plot (plot below shows only one PCA factor).
![Comp3](https://user-images.githubusercontent.com/88196987/232825772-94c5f075-c0b3-448e-a394-b05c01a221d2.jpeg)

I then select the associations that cross the significance threshold for multiple comparisons and plot those associations within the brain (see right panel for emergent associations)

![Comp2](https://user-images.githubusercontent.com/88196987/232826504-655f6a44-11eb-4ac2-bfd9-caec095617a3.JPG)
