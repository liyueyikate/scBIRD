#include <Rcpp.h>
#include "predict_header.h"


using namespace std;
//[[Rcpp::export]]
Rcpp::StringVector getLoci(Rcpp::CharacterVector libfile0="./src/human_hg19_model.bin"){
  std::string sep_str = ".";

  // initialize cpp data structure to hold rcpp data structure
  //char infile[255];
  char libfile[255];

  //copy rcpp data structure to cpp data strucutre
  //string infile_tmp = Rcpp::as<string>(infile0);
  //strcpy(infile, infile_tmp.c_str());
  string libfile_tmp = Rcpp::as<string>(libfile0);
  strcpy(libfile, libfile_tmp.c_str());

  Exondata exonin;
  //Exondata *indata = &exonin;
  int predictor_size,cluster_size, var_size, loci_size, bin_size, DH_num1, DH_num2, DH_num3;

  //read in the bin model file and get parameters for the regression model
  if (ReadPar(libfile, loci_size, predictor_size, cluster_size, bin_size, var_size, DH_num1, DH_num2, DH_num3))
  {
    return 1;
  }

  double *quantile_in = new double[predictor_size];
  int *cluster_idx = new int[predictor_size];
  double *exon_mean = new double[predictor_size];
  double *exon_sd = new double[predictor_size];
  double *DNase_mean = new double[loci_size];
  double *DNase_sd = new double[loci_size];
  double **coef = new double *[var_size];

  for (int j = 0; j < var_size; j++)
  {
    coef[j] = new double[loci_size]();
  }

  int **pre_idx = new int *[var_size];
  for (int j = 0; j < var_size; j++)
  {
    pre_idx[j] = new int[loci_size]();
  }

  char **select_loci = new char *[loci_size];
  for (int j = 0; j < loci_size; j++)
  {
    select_loci[j] = new char[30];
  }

  char **TC_id = new char *[predictor_size];
  for (int j = 0; j < predictor_size; j++)
  {
    TC_id[j] = new char[30];
  }

  //DH cluster model
  double **dis_matrix = new double *[3];
  for (int j = 0; j < 3; j++)
  {
    dis_matrix[j] = new double[loci_size]();
  }

  int **DH_cluster = new int *[3];
  for (int j = 0; j < 3; j++)
  {
    DH_cluster[j] = new int[loci_size]();
  }

  double **DH_coef1 = new double *[var_size];
  for (int j = 0; j < var_size; j++)
  {
    DH_coef1[j] = new double[DH_num1]();
  }

  double **DH_coef2 = new double *[var_size];
  for (int j = 0; j < var_size; j++)
  {
    DH_coef2[j] = new double[DH_num2]();
  }

  double **DH_coef3 = new double *[var_size];
  for (int j = 0; j < var_size; j++)
  {
    DH_coef3[j] = new double[DH_num3]();
  }

  int **DH_pre_idx1 = new int *[var_size];
  for (int j = 0; j < var_size; j++)
  {
    DH_pre_idx1[j] = new int[DH_num1]();
  }

  int **DH_pre_idx2 = new int *[var_size];
  for (int j = 0; j < var_size; j++)
  {
    DH_pre_idx2[j] = new int[DH_num2]();
  }

  int **DH_pre_idx3 = new int *[var_size];
  for (int j = 0; j < var_size; j++)
  {
    DH_pre_idx3[j] = new int[DH_num3]();
  }


  //read in model file
  if (ReadinModel(libfile, quantile_in, exon_mean, exon_sd, coef, DNase_mean, DNase_sd, pre_idx, TC_id, cluster_idx, select_loci, predictor_size, var_size, loci_size, dis_matrix, DH_cluster, DH_coef1, DH_coef2, DH_coef3, DH_pre_idx1, DH_pre_idx2, DH_pre_idx3, DH_num1, DH_num2, DH_num3))
  {
    return 1;
  }

  //Copy all loci to the all_loci_r vector and return
  Rcpp::StringVector all_loci_r(loci_size);

  for (int i = 0; i < loci_size; i++)
  {
    all_loci_r[i] = select_loci[i];
  }

  return(all_loci_r);

  //release memory.
  ReleaseExondata(exonin);
  delete[] quantile_in;
  delete[] cluster_idx;

  for (int i = 0; i < var_size; i++)
  {
    delete[] coef[i];
    delete[] pre_idx[i];
    delete[] DH_coef1[i];
    delete[] DH_coef2[i];
    delete[] DH_coef3[i];
    delete[] DH_pre_idx1[i];
    delete[] DH_pre_idx2[i];
    delete[] DH_pre_idx3[i];
  }
  for (int i = 0; i < 3; i++)
  {
    delete[] dis_matrix[i];
    delete[] DH_cluster[i];
  }
  for (int i = 0; i < predictor_size; i++)
  {
    delete[] TC_id[i];
  }
  for (int i = 0; i < loci_size; i++)
  {
    delete[] select_loci[i];
  }
  delete[] coef;
  delete[] pre_idx;
  delete[] DH_coef1;
  delete[] DH_coef2;
  delete[] DH_coef3;
  delete[] DH_pre_idx1;
  delete[] DH_pre_idx2;
  delete[] DH_pre_idx3;
  delete[] TC_id;
  delete[] select_loci;
  delete[] dis_matrix;
  delete[] DH_cluster;
}

//[[Rcpp::export]]
Rcpp::NumericMatrix bird_select_loci(Rcpp::CharacterVector infile0, int loci_size_tmp, Rcpp::NumericVector select_loci_tmp,Rcpp::CharacterVector libfile0="./src/human_hg19_model.bin",
                                     int locus_model=0, double up_bound=14, int match_mode=0,int write_flag = 0){

  std::string sep_str = ".";
  // initialize cpp data structure to hold rcpp data structure
  char infile[255];
  char libfile[255];

  //copy rcpp data structure to cpp data strucutre
  string infile_tmp = Rcpp::as<string>(infile0);
  strcpy(infile, infile_tmp.c_str());
  string libfile_tmp = Rcpp::as<string>(libfile0);
  strcpy(libfile, libfile_tmp.c_str());

  Exondata exonin;
  Exondata *indata = &exonin;
  int predictor_size, sample_size, cluster_size, var_size, loci_size, bin_size, DH_num1, DH_num2, DH_num3;

  //read in the bin model file and get parameters for the regression model
  if (ReadPar(libfile, loci_size, predictor_size, cluster_size, bin_size, var_size, DH_num1, DH_num2, DH_num3))
  {
    return 1;
  }

  double *quantile_in = new double[predictor_size];
  int *cluster_idx = new int[predictor_size];
  double *exon_mean = new double[predictor_size];
  double *exon_sd = new double[predictor_size];
  double *DNase_mean = new double[loci_size];
  double *DNase_mean_tmp = new double[loci_size_tmp];
  double *DNase_sd = new double[loci_size];
  double *DNase_sd_tmp = new double[loci_size_tmp];
  double **coef = new double *[var_size];
  double **coef_tmp = new double *[var_size];

  for (int j = 0; j < var_size; j++)
  {
    coef[j] = new double[loci_size]();
  }

  for (int j = 0; j < var_size; j++)
  {
    coef_tmp[j] = new double[loci_size_tmp]();
  }

  int **pre_idx = new int *[var_size];
  int **pre_idx_tmp = new int *[var_size];
  for (int j = 0; j < var_size; j++)
  {
    pre_idx[j] = new int[loci_size]();
  }

  for (int j = 0; j < var_size; j++)
  {
    pre_idx_tmp[j] = new int[loci_size_tmp]();
  }

  char **select_loci = new char *[loci_size];
  char **select_loci_tmp2 = new char *[loci_size_tmp]();
  for (int j = 0; j < loci_size; j++)
  {
    select_loci[j] = new char[30];
  }

  char **TC_id = new char *[predictor_size];
  for (int j = 0; j < predictor_size; j++)
  {
    TC_id[j] = new char[30];
  }

  //DH cluster model
  double **dis_matrix = new double *[3];
  double **dis_matrix_tmp = new double *[3];
  for (int j = 0; j < 3; j++)
  {
    dis_matrix[j] = new double[loci_size]();
  }

  for (int j = 0; j < 3; j++)
  {
    dis_matrix_tmp[j] = new double[loci_size_tmp]();
  }

  int **DH_cluster = new int *[3];
  int **DH_cluster_tmp = new int *[3];
  for (int j = 0; j < 3; j++)
  {
    DH_cluster[j] = new int[loci_size]();
  }

  for (int j = 0; j < 3; j++)
  {
    DH_cluster_tmp[j] = new int[loci_size_tmp]();
  }
  double **DH_coef1 = new double *[var_size];
  for (int j = 0; j < var_size; j++)
  {
    DH_coef1[j] = new double[DH_num1]();
  }

  double **DH_coef2 = new double *[var_size];
  for (int j = 0; j < var_size; j++)
  {
    DH_coef2[j] = new double[DH_num2]();
  }

  double **DH_coef3 = new double *[var_size];
  for (int j = 0; j < var_size; j++)
  {
    DH_coef3[j] = new double[DH_num3]();
  }

  int **DH_pre_idx1 = new int *[var_size];
  for (int j = 0; j < var_size; j++)
  {
    DH_pre_idx1[j] = new int[DH_num1]();
  }

  int **DH_pre_idx2 = new int *[var_size];
  for (int j = 0; j < var_size; j++)
  {
    DH_pre_idx2[j] = new int[DH_num2]();
  }

  int **DH_pre_idx3 = new int *[var_size];
  for (int j = 0; j < var_size; j++)
  {
    DH_pre_idx3[j] = new int[DH_num3]();
  }

  //read in model file
  if (ReadinModel(libfile, quantile_in, exon_mean, exon_sd, coef, DNase_mean, DNase_sd, pre_idx, TC_id, cluster_idx, select_loci, predictor_size, var_size, loci_size, dis_matrix, DH_cluster, DH_coef1, DH_coef2, DH_coef3, DH_pre_idx1, DH_pre_idx2, DH_pre_idx3, DH_num1, DH_num2, DH_num3))
  {
    return 1;
  }

  //create tmp param only for the selected loci
  for (int k=0; k< loci_size;k++){
    for (int j=0; j<loci_size_tmp;j++) {
      if (k==select_loci_tmp[j]){
        DNase_mean_tmp[j]=DNase_mean[k];
      }
    }
  }

  for (int k=0; k< loci_size;k++){
    for (int j=0; j<loci_size_tmp;j++) {
      if (k==select_loci_tmp[j]){
        select_loci_tmp2[j]=select_loci[k];
      }
    }
  }

  for (int k=0; k< loci_size;k++){
    for (int j=0; j<loci_size_tmp;j++){
      if (k==select_loci_tmp[j]){
        DNase_sd_tmp[j]= DNase_sd[k];
      }
    }
  }

  for (int k=0; k<var_size;k++){
    for (int j=0; j<loci_size_tmp;j++){
      for (int m=0;m<loci_size;m++){
        if (m==select_loci_tmp[j]){
          coef_tmp[k][j]=coef[k][m];
        }
      }
    }
  }

  for (int k=0; k<var_size;k++){
    for (int j=0; j<loci_size_tmp;j++){
      for (int m=0;m<loci_size;m++){
        if (m==select_loci_tmp[j]){
          pre_idx_tmp[k][j]= pre_idx[k][m];
        }
      }
    }
  }

  for (int k=0; k<3;k++){
    for (int j=0; j<loci_size_tmp;j++){
      for (int m=0;m<loci_size;m++){
        if (m==select_loci_tmp[j]){
          dis_matrix_tmp[k][j]= dis_matrix[k][m];
        }
      }
    }
  }

  for (int k=0; k<3;k++){
    for (int j=0; j<loci_size_tmp;j++){
      for (int m=0;m<loci_size;m++){
        if (m==select_loci_tmp[j]){
          DH_cluster_tmp[k][j]= DH_cluster[k][m];
        }
      }
    }
  }

  DNase_mean=DNase_mean_tmp;
  DNase_sd=DNase_sd_tmp;
  coef=coef_tmp;
  pre_idx=pre_idx_tmp;
  select_loci=select_loci_tmp2;
  dis_matrix=dis_matrix_tmp;
  DH_cluster=DH_cluster_tmp;
  loci_size=loci_size_tmp;

  if (ReadinExon(infile, indata))
  {
    return 1;
  }

  std::cout <<indata << std::endl;
  int *match_idx = new int[predictor_size];

  //match and check gene expression data
  std::cout << "Matching the input gene expression data." << std::endl;
  if(match_mode == 1)
  {
    //use exact match
    if (MatchExon(TC_id, indata, predictor_size, match_idx))
    {
      std::cout << "Gene expression data format incorrect: No gene id matches the library file." << std::endl;
      std::cout << "Please check sample input file for reference." << std::endl;
      return 1;
    }
  }
  else
  {
    //match id before separate character
    if (MatchExon_sep(TC_id, indata, predictor_size, match_idx, sep_str))
    {
      std::cout << "Gene expression data format incorrect: No gene id matches the library file." << std::endl;
      std::cout << "Please check sample input file for reference." << std::endl;
      return 1;
    }
  }


  std::cout << "Processing data..." << std::endl;
  sample_size = (int)exonin.sample_name.size();
  double ** data_norm = new double *[sample_size];
  double ** data_mean = new double *[sample_size];
  double ** output = new double *[sample_size];
  double ** DH_pre1 = new double *[sample_size];
  double ** DH_pre2 = new double *[sample_size];
  double ** DH_pre3 = new double *[sample_size];

  for (int i = 0; i < sample_size; i++)
  {
    data_norm[i] = new double[predictor_size]();
    data_mean[i] = new double[cluster_size]();
    output[i] = new double[loci_size]();
    DH_pre1[i] = new double[DH_num1]();
    DH_pre2[i] = new double[DH_num2]();
    DH_pre3[i] = new double[DH_num3]();
  }

  for (int i = 0; i < sample_size; i++)
  {
    for (int j = 0; j < predictor_size; j++)
    {
      if(match_idx[j] != -1)
      {
        data_norm[i][j] = log2(exonin.data[match_idx[j]][i] + 1);
      }
      else
      {
        data_norm[i][j] = 0;
      }
    }
  }

  for (int i = 0; i < sample_size; i++)
  {
    QuantileNorm(data_norm[i], quantile_in, predictor_size);
  }

  StandardizeRow(data_norm, exon_mean, exon_sd, predictor_size, sample_size); //standardize gene expression data
  ClusterMean(data_norm, data_mean, cluster_idx, predictor_size, cluster_size, sample_size); //get gene expression cluster mean;

  if(locus_model != 1)
  {
    //Locus level prediction
    Regression(data_mean, output, coef, pre_idx, var_size, loci_size, sample_size);

    //DH cluster level prediction
    Regression(data_mean, DH_pre1, DH_coef1, DH_pre_idx1, var_size, DH_num1, sample_size);
    Regression(data_mean, DH_pre2, DH_coef2, DH_pre_idx2, var_size, DH_num2, sample_size);
    Regression(data_mean, DH_pre3, DH_coef3, DH_pre_idx3, var_size, DH_num3, sample_size);

    //Combine result from locus level and DH cluster level
    ModelAverage(output, DH_pre1, DH_pre2, DH_pre3, dis_matrix, DH_cluster, loci_size, sample_size);

    std::cout << "Using full model for prediction." << std::endl;

  }
  else
  {
    //Locus level prediction
    Regression(data_mean, output, coef, pre_idx, var_size, loci_size, sample_size);
    std::cout << "Using locus-level model for prediction." << std::endl;

  }

  //convert predicted value back to original scale
  StandardizeRow_r(output, DNase_mean, DNase_sd, loci_size, sample_size);

  // return output;
  Rcpp::NumericMatrix output_r(sample_size,loci_size);

  for (int i = 0; i < loci_size; i++)
  {
    for (int j = 0; j < sample_size; j++)
    {
      output_r(j,i) = output[j][i];
    }
  }


  return transpose(output_r);

  //release memory.
  ReleaseExondata(exonin);
  delete[] quantile_in;
  delete[] cluster_idx;
  for (int i = 0; i < sample_size; i++)
  {
    delete[] data_norm[i];
    delete[] data_mean[i];
    delete[] output[i];
    delete[] DH_pre1[i];
    delete[] DH_pre2[i];
    delete[] DH_pre3[i];
  }
  for (int i = 0; i < var_size; i++)
  {
    delete[] coef[i];
    delete[] coef_tmp[i];
    delete[] pre_idx[i];
    delete[] DH_coef1[i];
    delete[] DH_coef2[i];
    delete[] DH_coef3[i];
    delete[] DH_pre_idx1[i];
    delete[] DH_pre_idx2[i];
    delete[] DH_pre_idx3[i];
  }
  for (int i = 0; i < 3; i++)
  {
    delete[] dis_matrix[i];
    delete[] dis_matrix_tmp[i];
    delete[] DH_cluster[i];
    delete[] DH_cluster_tmp[i];
  }

  for (int i = 0; i < predictor_size; i++)
  {
    delete[] TC_id[i];
  }
  for (int i = 0; i < loci_size; i++)
  {
    delete[] select_loci[i];
  }

  for (int i = 0; i < loci_size_tmp; i++)
  {
    delete[] select_loci_tmp2[i];
  }
  delete[] data_norm;
  delete[] data_mean;
  delete[] output;
  delete[] DH_pre1;
  delete[] DH_pre2;
  delete[] DH_pre3;
  delete[] coef;
  delete[] coef_tmp;
  delete[] pre_idx;
  delete[] DH_coef1;
  delete[] DH_coef2;
  delete[] DH_coef3;
  delete[] DH_pre_idx1;
  delete[] DH_pre_idx2;
  delete[] DH_pre_idx3;
  delete[] TC_id;
  delete[] select_loci;
  delete[] select_loci_tmp2;
  delete[] dis_matrix;
  delete[] dis_matrix_tmp;
  delete[] DH_cluster;
  delete[] DH_cluster_tmp;
}




