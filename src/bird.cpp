#include <Rcpp.h>
#include "predict_header.h"


using namespace std;



// [[Rcpp::export]]

Rcpp::NumericMatrix bird(Rcpp::CharacterVector infile0, Rcpp::CharacterVector outfile0, Rcpp::CharacterVector outfile_expr0, Rcpp::CharacterVector libfile0="./src/human_hg19_model.bin",
 int locus_model=0, double up_bound=14, int match_mode=0,int write_flag = 0){

  std::string sep_str = ".";

  // initialize cpp data structure to hold rcpp data structure
  char infile[255];
  char outfile[255];
  char outfile_expr[255];
  char libfile[255];

  //copy rcpp data structure to cpp data strucutre
  string infile_tmp = Rcpp::as<string>(infile0);
  strcpy(infile, infile_tmp.c_str());
  string outfile_tmp = Rcpp::as<string>(outfile0);
  strcpy(outfile, outfile_tmp.c_str());
  string outfile_expr_tmp = Rcpp::as<string>(outfile_expr0);
  strcpy(outfile_expr, outfile_expr_tmp.c_str());
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


  //check chromosome
 // char **select_loci0 = select_loci;
  //std::cout <<select_loci0[2] << std::endl;
  //std::cout <<loci_size << std::endl;


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

  WriteExpr(data_norm, exonin.sample_name, TC_id, outfile_expr, predictor_size, sample_size);

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

    /*
    std::cout <<DH_pre1[1][1] << std::endl;
    std::cout <<DH_pre2[1][1] << std::endl;
    std::cout <<DH_pre3[1][1] << std::endl;
    std::cout <<output[1][1] << std::endl;
     */
  }
  else
  {
    //Locus level prediction
    Regression(data_mean, output, coef, pre_idx, var_size, loci_size, sample_size);
    std::cout << "Using locus-level model for prediction." << std::endl;

  }

  //convert predicted value back to original scale
  StandardizeRow_r(output, DNase_mean, DNase_sd, loci_size, sample_size);

  //write output file
  std::cout << "Writing output file..." << std::endl;
  if (WriteWIG(output, select_loci, exonin.sample_name, outfile, bin_size, loci_size, sample_size, write_flag, up_bound))
  {
    return 1;
  }


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
  delete[] data_norm;
  delete[] data_mean;
  delete[] output;
  delete[] DH_pre1;
  delete[] DH_pre2;
  delete[] DH_pre3;
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
  /*
  return 0;
   */
}
