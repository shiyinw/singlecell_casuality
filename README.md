# Single-Cell + PPI network -> Unsupervised Casuality Discovery
Course Project for "Hot Topics in Computational Biology" in 2019 Autumn

## Advice by TA: Infer gene regulatory network by using graphical

通过各种模型来构建网络

Condition-adaptive fused graphical lasso (CFGL): An adaptive procedure for inferring condition-specific gene co-expression network

Identifying gene regulatory network rewiring using latent differential graphical models

## Related Works

#### [immunology2017] Single-Cell Genomics: Approaches and Utility in Immunology

https://www.cell.com/trends/immunology/fulltext/S1471-4906(16)30205-8

综述和工具介绍，学习术语

#### [best2019] Current best practices in single‐cell RNA‐seq analysis: a tutorial

https://www.embopress.org/doi/full/10.15252/msb.20188746

综述类文章

提供了可以上手的实例代码：https://github.com/shiyinw/single-cell-tutorial

#### [integration2019] Efficient integration of heterogeneous single-cell transcriptomes using Scanorama

https://www.nature.com/articles/s41587-019-0113-3#data-availability

Bonnie Burger的数据集聚合成果，也有很多的数据集，会比较新

#### [integration2018] Integrating single-cell transcriptomic data across different conditions, technologies, and species

https://www.nature.com/articles/nbt.4096?draft=journal#Sec12

#### [cluster2018] Single cell clustering based on cell-pair differentiability correlation and variance analysis

https://academic.oup.com/bioinformatics/article/34/21/3684/4996592?searchresult=1

做clustering的，学习如何给cell进行概率建模

#### [gen2018] Deep generative modeling for single-cell transcriptomics (Nature Methods, Nov 2018)

https://www.nature.com/articles/s41592-018-0229-2

数据部分特别全，从中选择带有gene expression的数据集

#### [ode2017] SCODE: an efficient regulatory network inference algorithm from single-cell RNA-Seq during differentiation

https://academic.oup.com/bioinformatics/article/33/15/2314/3100331?searchresult=1#118776066

一篇用ODE来分析single-cell RNA-Seq来生成网络的工作，方法上有借鉴意义

#### [model2019] Data Integration of Hybrid Microarray and Single Cell Expression Data to Enhance Gene Network Inference

http://www.eurekaselect.com/node/168772/viewhtml/U2s5mVmhpVaz0VkQrlRGoRU1b2UTmBrVrk9BkTDBxOdkxNTVgWg4cdk1ltVF8F2ceU15COWh8vZECdjMEXNMpazBqOQ2aVTVxTh36dE10UcVaF0TdWR5qMHlyoTl6RVahnVhtSHddSdGzJqS9Hg4fYmRBHVnjM0ZqENwtOW9mkbEccxca3JmcSHdo5MVzpUZE0V3fWWFbqUXrpkTajJVrQnlsNMmctVegU5s8ekZcrdFjptV3jVrgWkQuxbGApZb1zJJyd2ptaVGsttdi01UAZ0JgsWmdZue9Dh0tWVhfKaToBhVb3BOuc1pEkUVFRjVellUFY1ZrZNXzx8VrldspeGx3SV2cxnd2VpFjQ1V70UXoVWTqkU0lYzJ4SZlapaao2JW3SlFDmVEw15YcVJEzTW41xUmBpjYyzFOleURjWMkzVMVkFFj4M01e6M0ol0TiXBqbTTFtuTmepRdqHpNfajlgORVwJ4MplZFsZGVpBVGkNWWiTNUdY1ZlZ

如何给网络进行建模，拓扑结构->概率

## Code

https://github.com/shiyinw/single-cell-tutorial



