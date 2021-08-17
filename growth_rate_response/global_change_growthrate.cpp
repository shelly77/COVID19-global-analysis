
//The only difference between this file and "stateSpaceGlobalAutocor.cpp"
//is that this file removes strigency index from the model, since it is correlated
//with population mobility


#include <TMB.hpp>

// template <class Type>
//   Type f(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);} //convert x to range (0,1)
// 
// 
// // template<class Type> // define generalised poisson distribution
// // inline Type dgpois(const Type &x, const Type &lambda, const Type &omega, int give_log=0)
// // {
//   //   Type lambdastar = (Type(1) - omega)*lambda + omega*x;
//   //   Type logres = log(1 - omega) + log(lambda) + (x - 1)*log(lambdastar) -
//     //     lgamma(x + 1) - lambdastar;
//   //   if (give_log) return logres; else return exp(logres);
//   // }
// 
// 
// template<class Type> // define zero-inflated binomial distribution
// inline Type dzibinom(const Type &x, const Type &size, const Type &p, const Type & zip,
//                      int give_log=0)
// {
//   Type logres;
//   if (x==Type(0)) logres=log(zip + (Type(1)-zip)*dbinom(x, size, p, false)); 
//   else logres=log(Type(1)-zip) + dbinom(x, size, p, true);
//   if (give_log) return logres; else return exp(logres);
// }
// 
template<class Type> // define beta-binomial distribution
inline Type dbetabinom(const Type &x, const Type &size, const Type &alpha, const Type &beta,
                       int give_log=0) // size is the total number of trials and x is the number of success
{
  Type logres =
    lfactorial(size)
  - lfactorial(x)
  - lfactorial(size - x)
  + lgamma(alpha + x)
  + lgamma(alpha + beta)
  - lgamma(alpha + beta + size)
  - lgamma(alpha)
  - lgamma(beta)
  + lgamma(beta + size - x+ 1e-100);
  if (give_log) return logres; else return exp(logres);
}
// 
// template<class Type> // define zero-inflated beta-binomial distribution
// inline Type dzibetabinom(const Type &x, const Type &size, const Type &alpha, const Type &beta, const Type & zip,
//                          int give_log=0) // size is the total number of trials and x is the number of success
// {
//   Type logres;
//   if (x==Type(0)) logres=log(zip + (Type(1)-zip)*dbetabinom(x, size, alpha, beta, false)); 
//   else logres=log(Type(1)-zip) + dbetabinom(x, size, alpha, beta, true);
//   if (give_log) return logres; else return exp(logres);
// }


template<class Type>
Type objective_function<Type>::operator() (){
  
  using namespace density;
  
  
  DATA_VECTOR(grorate_change); 
  
  DATA_VECTOR(temp);
  DATA_VECTOR(temp_sq);
  DATA_VECTOR(temp_lag1);
  DATA_VECTOR(temp_lag2);
  DATA_VECTOR(windsp);
  DATA_VECTOR(windsp_lag1);
  DATA_VECTOR(windsp_lag2);
  DATA_VECTOR(humidity);
  DATA_VECTOR(humidity_lag1);
  DATA_VECTOR(humidity_lag2);
  DATA_VECTOR(uv);
  DATA_VECTOR(uv_lag1);
  DATA_VECTOR(uv_lag2);
  
  // 
  DATA_VECTOR(population);
  DATA_VECTOR(dens);
  DATA_VECTOR(mage);
  // 
  DATA_VECTOR(gdp);
  // 
  DATA_VECTOR(newtests);
  // 
  
  DATA_VECTOR(mobility);
  DATA_VECTOR(mobility_lag1);
  DATA_VECTOR(mobility_lag2);
  
  DATA_VECTOR(debtrelief);
  DATA_VECTOR(debtrelief_lag1);
  DATA_VECTOR(debtrelief_lag2);
  
  DATA_VECTOR(healthinvestment);
  DATA_VECTOR(healthinvestment_lag1);
  DATA_VECTOR(healthinvestment_lag2);
  
  DATA_VECTOR(contatracing);
  DATA_VECTOR(contatracing_lag1);
  DATA_VECTOR(contatracing_lag2);
  
  DATA_VECTOR(cases_days_since_first);
  DATA_VECTOR(cases_days_since_first_sq);
  
  DATA_FACTOR(spatioLevelfactor);
  
  DATA_VECTOR(index);
  
  PARAMETER(beta0); 
  
  PARAMETER(beta_temp);
  PARAMETER(beta_temp_sq);
  PARAMETER(beta_temp_lag1);
  PARAMETER(beta_temp_lag2);
  
  PARAMETER(beta_windsp);
  PARAMETER(beta_windsp_lag1);
  PARAMETER(beta_windsp_lag2);
  
  PARAMETER(beta_humidity);
  PARAMETER(beta_humidity_lag1);
  PARAMETER(beta_humidity_lag2);
  
  PARAMETER(beta_uv);
  PARAMETER(beta_uv_lag1);
  PARAMETER(beta_uv_lag2);
  // // 
  PARAMETER(beta_population);
  PARAMETER(beta_dens);
  PARAMETER(beta_mage);
  // 
  PARAMETER(beta_gdp);
  // 
  PARAMETER(beta_newtests);
  
  PARAMETER(beta_mobility);
  PARAMETER(beta_mobility_lag1);
  PARAMETER(beta_mobility_lag2);
  
  PARAMETER(beta_debtrelief);
  PARAMETER(beta_debtrelief_lag1);
  PARAMETER(beta_debtrelief_lag2);
  
  PARAMETER(beta_healthinvestment);
  PARAMETER(beta_healthinvestment_lag1);
  PARAMETER(beta_healthinvestment_lag2);
  
  PARAMETER(beta_contatracing);
  PARAMETER(beta_contatracing_lag1);
  PARAMETER(beta_contatracing_lag2);
  
  PARAMETER(beta_days);
  PARAMETER(beta_days_sq);
  
  //random spatioLevel effects
  PARAMETER_VECTOR(mu_spatioLevel);
  PARAMETER(logbeta_spatioLevel);   Type beta_spatioLevel = exp(logbeta_spatioLevel);  ADREPORT(exp(2*beta_spatioLevel));
  
  //random slopes of days
  PARAMETER_VECTOR(mu_spatioLevel_days);
  PARAMETER(logsd_randm_slop_days);   Type sd_randm_slop_days = exp(logsd_randm_slop_days);  ADREPORT(exp(2*logsd_randm_slop_days));
  
  //random slopes of squared days
  PARAMETER_VECTOR(mu_spatioLevel_sqdays);
  PARAMETER(logsd_randm_slop_sqdays);   Type sd_randm_slop_sqdays = exp(logsd_randm_slop_sqdays); ADREPORT(exp(2*logsd_randm_slop_sqdays));
  
  //random slopes on mobility
  PARAMETER_VECTOR(mu_spatioLevel_mobility); 
  PARAMETER(logsd_randm_slop_mobility);   Type sd_randm_slop_mobility= exp(logsd_randm_slop_mobility); ADREPORT(exp(2*logsd_randm_slop_mobility));
  
  //random slopes of one week lagged mobility
  PARAMETER_VECTOR(mu_spatioLevel_mobility_lag1);
  PARAMETER(logsd_randm_slop_mobility_lag1);   Type sd_randm_slop_mobility_lag1 = exp(logsd_randm_slop_mobility_lag1); ADREPORT(exp(2*logsd_randm_slop_mobility_lag1));
  
  //random slopes of two weeks lagged mobility
  PARAMETER_VECTOR(mu_spatioLevel_mobility_lag2); 
  PARAMETER(logsd_randm_slop_mobility_lag2);   Type sd_randm_slop_mobility_lag2 = exp(logsd_randm_slop_mobility_lag2); ADREPORT(exp(2*logsd_randm_slop_mobility_lag2));
  
  //random slopes on temperature
  PARAMETER_VECTOR(mu_spatioLevel_temp);
  PARAMETER(logsd_randm_slop_temp);   Type sd_randm_slop_temp = exp(logsd_randm_slop_temp); ADREPORT(exp(2*logsd_randm_slop_temp));
  
  //random slopes on squared temperature
  PARAMETER_VECTOR(mu_spatioLevel_sqtemp);
  PARAMETER(logsd_randm_slop_sqtemp);   Type sd_randm_slop_sqtemp = exp(logsd_randm_slop_sqtemp); ADREPORT(exp(2*logsd_randm_slop_temp));
  
  // parameters in ARCH(2)
  PARAMETER(logalpha0); Type alpha0 = exp(logalpha0); ADREPORT(exp(2*logalpha0));
  PARAMETER(logalpha1); Type alpha1 = exp(logalpha1); ADREPORT(exp(2*logalpha1));
  PARAMETER(logalpha2); Type alpha2 = exp(logalpha2); ADREPORT(exp(2*logalpha2));
  
  //standard deviation of autoregressive error
  //PARAMETER(logsd_error_measure);  Type sd_error_measure = exp(logsd_error_measure); ADREPORT(sd_error_measure);
  
  
  int nObs  = grorate_change.size();
  //int nspatioLevel = NLEVELS(spatioLevel);
  
  Type nll=Type(0);
  
  
  //vector<Type> error_sate(nObs);
  
  //nll -= dnorm(,Type(0),Type(1),1));
  
  // for (int i=1; i<nObs; i++){
  
  //     for (int t=0; t<nYear; t++){
  //       mu(t,0)= mu_alpha(t);
  //       mu(t,1)= mu_theta(t);
  //       mu(t,2)= mu_omega(t);
  //     }
  //     
  //     // Rescaling such such that mu_alpha, mu_theta and mu_omega are standard normal  
  //     vector<Type> sds(3);
  //     sds(0) = 1/sqrt(Gamma0(0,0));
  //     sds(1) = 1/sqrt(Gamma0(1,1));
  //     sds(2) = 1/sqrt(Gamma0(2,2));
  //     
  //     nll += VECSCALE(MVNORM(Gamma0),sds)(mu.row(0));
  //     
  // 
  //contribution of random spatioLevel effects to likelihood
  if (CppAD::Variable(mu_spatioLevel(0)))
    nll-= sum(dnorm(mu_spatioLevel,Type(0),Type(1),1)); //iid N(0,1) random year effect
  
  if (CppAD::Variable(mu_spatioLevel_days(0)))
    nll-= sum(dnorm(mu_spatioLevel_days,Type(0),sd_randm_slop_days,1)); //iid N(0,1) random year effect
  
  if (CppAD::Variable(mu_spatioLevel_sqdays(0)))
    nll-= sum(dnorm(mu_spatioLevel_sqdays,Type(0),sd_randm_slop_sqdays,1)); //iid N(0,1) random year effect
  
  if (CppAD::Variable(mu_spatioLevel_mobility(0)))
    nll-= sum(dnorm(mu_spatioLevel_mobility,Type(0),sd_randm_slop_mobility,1)); //iid N(0,1) random year effect
  
  if (CppAD::Variable(mu_spatioLevel_mobility_lag1(0)))
    nll-= sum(dnorm(mu_spatioLevel_mobility_lag1,Type(0),sd_randm_slop_mobility_lag1,1)); //iid N(0,1) random year effect
  
  if (CppAD::Variable(mu_spatioLevel_mobility_lag2(0)))
    nll-= sum(dnorm(mu_spatioLevel_mobility_lag2,Type(0),sd_randm_slop_mobility_lag2,1)); //iid N(0,1) random year effect
  
  if (CppAD::Variable(mu_spatioLevel_temp(0)))
    nll-= sum(dnorm(mu_spatioLevel_temp,Type(0),sd_randm_slop_temp,1)); //iid N(0,1) random year effect
  
  if (CppAD::Variable(mu_spatioLevel_sqtemp(0)))
    nll-= sum(dnorm(mu_spatioLevel_sqtemp,Type(0),sd_randm_slop_sqtemp,1)); //iid N(0,1) random year effect
  
  
  vector<Type> state(nObs);
  vector<Type> error(nObs);
  
  for (int i=0; i<nObs; i++){
    
    // Type etaEnv = beta0;  
    Type eta = beta0;
    
    if (CppAD::Variable(beta_temp))
      eta += beta_temp*temp(i);
    if (CppAD::Variable(beta_temp_sq))
      eta += beta_temp_sq*temp_sq(i);
    if (CppAD::Variable(beta_temp_lag1))
      eta += beta_temp_lag1*temp_lag1(i);
    if (CppAD::Variable(beta_temp_lag2))
      eta += beta_temp_lag2*temp_lag2(i);
    
    if (CppAD::Variable(beta_windsp))
      eta += beta_windsp*windsp(i);
    if (CppAD::Variable(beta_windsp_lag1))
      eta += beta_windsp_lag1*windsp_lag1(i);
    if (CppAD::Variable(beta_windsp_lag2))
      eta += beta_windsp_lag2*windsp_lag2(i);
    
    if (CppAD::Variable(beta_humidity))
      eta += beta_humidity*humidity(i);
    if (CppAD::Variable(beta_humidity_lag1))
      eta += beta_humidity_lag1*humidity_lag1(i);
    if (CppAD::Variable(beta_humidity_lag2))
      eta += beta_humidity_lag2*humidity_lag2(i);
    
    if (CppAD::Variable(beta_uv))
      eta += beta_uv*uv(i);
    if (CppAD::Variable(beta_uv_lag1))
      eta += beta_uv_lag1*uv_lag1(i);
    if (CppAD::Variable(beta_uv_lag2))
      eta += beta_uv_lag2*uv_lag2(i);
    // 
    if (CppAD::Variable(beta_population))
      eta += beta_population*population(i);
    
    if (CppAD::Variable(beta_dens))
      eta += beta_dens*dens(i);
    // 
    if (CppAD::Variable(beta_mage))
      eta += beta_mage*mage(i);
    // 
    if (CppAD::Variable(beta_gdp))
      eta += beta_gdp*gdp(i);
    
    if (CppAD::Variable(beta_newtests))
      eta += beta_newtests*newtests(i);
    // 
    if (CppAD::Variable(beta_mobility))
      eta += beta_mobility*mobility(i);
    if (CppAD::Variable(beta_mobility_lag1))
      eta += beta_mobility_lag1*mobility_lag1(i);
    if (CppAD::Variable(beta_mobility_lag2))
      eta += beta_mobility_lag2*mobility_lag2(i);
    
    if (CppAD::Variable(beta_mobility))
      eta += beta_mobility*mobility(i);
    if (CppAD::Variable(beta_mobility_lag1))
      eta += beta_mobility_lag1*mobility_lag1(i);
    if (CppAD::Variable(beta_mobility_lag2))
      eta += beta_mobility_lag2*mobility_lag2(i);
    
    if (CppAD::Variable(beta_debtrelief))
      eta += beta_debtrelief*debtrelief(i);
    if (CppAD::Variable(beta_debtrelief_lag1))
      eta += beta_debtrelief_lag1*debtrelief_lag1(i);
    if (CppAD::Variable(beta_debtrelief_lag2))
      eta += beta_debtrelief_lag2*debtrelief_lag2(i);
    
    if (CppAD::Variable(beta_healthinvestment))
      eta += beta_healthinvestment*healthinvestment(i);
    if (CppAD::Variable(beta_healthinvestment_lag1))
      eta += beta_healthinvestment_lag1*healthinvestment_lag1(i);
    if (CppAD::Variable(beta_healthinvestment_lag2))
      eta += beta_healthinvestment_lag2*healthinvestment_lag2(i);
    
    if (CppAD::Variable(beta_contatracing))
      eta += beta_contatracing*contatracing(i);
    if (CppAD::Variable(beta_contatracing_lag1))
      eta += beta_contatracing_lag1*contatracing_lag1(i);
    if (CppAD::Variable(beta_contatracing_lag2))
      eta += beta_contatracing_lag2*contatracing_lag2(i);
    
    if (CppAD::Variable(beta_days))
      eta += beta_days*cases_days_since_first(i); 
    
    if (CppAD::Variable(beta_days_sq))
      eta +=beta_days_sq*cases_days_since_first_sq(i); 
    
    if (CppAD::Variable(mu_spatioLevel(0)))
      eta += beta_spatioLevel*mu_spatioLevel(spatioLevelfactor(i));
    
    if (CppAD::Variable(mu_spatioLevel_days(0)))
      eta += mu_spatioLevel_days(spatioLevelfactor(i))*cases_days_since_first(i);
    
    if (CppAD::Variable(mu_spatioLevel_sqdays(0)))
      eta += mu_spatioLevel_sqdays(spatioLevelfactor(i))*cases_days_since_first_sq(i);
    
    if (CppAD::Variable(mu_spatioLevel_mobility(0)))
      eta += mu_spatioLevel_mobility(spatioLevelfactor(i))*mobility(i);
    
    if (CppAD::Variable(mu_spatioLevel_mobility_lag1(0)))
      eta += mu_spatioLevel_mobility_lag1(spatioLevelfactor(i))*mobility_lag1(i);
    
    if (CppAD::Variable(mu_spatioLevel_mobility_lag2(0)))
      eta += mu_spatioLevel_mobility_lag2(spatioLevelfactor(i))*mobility_lag2(i);
    
    if (CppAD::Variable(mu_spatioLevel_temp(0)))
      eta += mu_spatioLevel_temp(spatioLevelfactor(i))*temp(i);
    
    if (CppAD::Variable(mu_spatioLevel_sqtemp(0)))
      eta += mu_spatioLevel_sqtemp(spatioLevelfactor(i))*temp_sq(i);
    
    
    state(i) = eta;
    
    error(i) = grorate_change(i) - state(i);
    
    // error_measure(i) = (newcasessmooth(i) - Gamma*state)/sd_error_measure;
    
    // error_measure(i) = (newcasessmooth(i) - state)/sd_error_measure;
    
    
    //  nll -= dnorm(error_measure(i), Type(0), Type(1), true);
    
    //nll -= dnorm(grorate_change(i), state(i), sd_error_measure, true);
    
  }
  
  
  for (int i=0; i<nObs; i++) {
    if(index(i) == 1) nll -= dnorm(error(i), Type(0), alpha0, true);
    else{
      Type sd = alpha0;
      if (CppAD::Variable(logalpha1) &  !CppAD::Variable(logalpha2))    sd = sqrt(alpha0+alpha1*pow(error(i-1),2));
      else if (CppAD::Variable(logalpha1) & CppAD::Variable(logalpha2))  sd = sqrt(alpha0+alpha1*pow(error(i-1),2)+alpha2*pow(error(i-2),2));
      nll -= dnorm(error(i), Type(0), sd, true);
    }
  }
  
  ADREPORT(state);
  ADREPORT(error);
  
  return nll;
}





// Type etaSoceco = etaSoceco
//   
//   if (CppAD::Variable(logtheta_year(s)))
//     etaEnv += theta_year(s)*mu_theta(t(i)); 
//   if (CppAD::Variable(theta_peakd(s)))
//     etaEnv +=  theta_peakd(s)*peakd(t(i));
//   if (CppAD::Variable(theta_clusize(s)))
//     etaEnv += theta_clusize(s)*y(i,0);  
//   if (CppAD::Variable(logtheta_year(s)))
//     etaEnv += theta_year(s)*mu_theta(t(i)); 



// Type mean = etaEnv + etaDemo + etaSoceco

//   nll -= dbinom(y(i,4),Type(1), survival(i,3), true);}









//   PARAMETER_VECTOR(logalpha_year); vector<Type> alpha_year = exp(logalpha_year); ADREPORT(alpha_year);
//   
//   
//   
//   PARAMETER_VECTOR(alpha_BCI2);
//   PARAMETER_VECTOR(alpha_BCI3);
//   PARAMETER_VECTOR(alpha_clusize);
//   PARAMETER_VECTOR(alpha_dens);
//   PARAMETER_VECTOR(logalpha_mother);  vector<Type> alpha_mother = exp(logalpha_mother); ADREPORT(alpha_mother);
//   PARAMETER_VECTOR(logalpha_year); vector<Type> alpha_year = exp(logalpha_year); ADREPORT(alpha_year);
//   
//   //when the cluth size is under directional selection, we use a glm with random intercepts and slopes 
//   //to model clutch size against laying date
//   // PARAMETER(alpha1_laydate); //fixed slope of directional selection for episode 1
//   // PARAMETER_VECTOR(alpha0_t); //random intercepts
//   // PARAMETER_VECTOR(alpha1_t); //random slopes
//   // //sdVecRanP is vector of standard deviations of random intercepts and slopes for episode 1
//   // PARAMETER_VECTOR(logSdRan);    vector<Type> sdVecRan = exp(logSdRan); ADREPORT(sdVecRan);
//   // //rhoRan is correlation between random intercepts and slopes for component P
//   // PARAMETER(rhoRan);  
//   
//   
//   PARAMETER_VECTOR(theta_0);
//   PARAMETER_VECTOR(theta_temp);
//   PARAMETER_VECTOR(theta_dens);
//   PARAMETER_VECTOR(theta_peakd);
//   PARAMETER_VECTOR(theta_clusize);
//   PARAMETER_VECTOR(logtheta_year); vector<Type> theta_year = exp(logtheta_year); ADREPORT(theta_year);
//   
//   PARAMETER_VECTOR(omega_0);
//   PARAMETER_VECTOR(omega_temp);
//   PARAMETER_VECTOR(logomega_year); vector<Type> omega_year = exp(logomega_year); ADREPORT(omega_year);
//   
//   PARAMETER_VECTOR(zip_0);
//   PARAMETER_VECTOR(beta_RLD);
//   PARAMETER_VECTOR(beta_dens);
//   PARAMETER_VECTOR(beta_diffSurvFem);
//   PARAMETER_VECTOR(beta_clusize);
//   PARAMETER_VECTOR(beta_peakd);
//   PARAMETER_VECTOR(beta_LayDate);
//   PARAMETER_VECTOR(logbeta_year);    vector<Type> beta_year = exp(logbeta_year); ADREPORT(beta_year);
//   
//   //random mother effects
//   PARAMETER_VECTOR(epsi_surv1);
//   PARAMETER_VECTOR(epsi_surv2);
//   PARAMETER_VECTOR(epsi_surv3);
//   PARAMETER_VECTOR(epsi_mothersurv);
//   //cormother are the correlations between the random mother effects
//   
//   //random year effects
//   PARAMETER_VECTOR(mu_alpha);
//   PARAMETER_VECTOR(mu_theta);
//   PARAMETER_VECTOR(mu_omega);
//   PARAMETER_VECTOR(mu_zip);
//   
//   //parCor are the correlations between the errors of mu_alpha,mu_theta and mu_omega
//   PARAMETER_VECTOR(parCor);
//   //Phi is the autocorrelation matrix of mu_alpha,mu_theta and mu_omega
//   PARAMETER_MATRIX(phi);  matrix<Type>  Phi =2/(1+exp(-2*phi.array()))-1; ADREPORT(Phi);
//   
//   //dispersion parameters for the zero-inflated negative bionomial
//   PARAMETER_VECTOR(log_disPar); vector<Type> disPar = 1 + 1e-10 + exp(log_disPar); ADREPORT(disPar);
//   
//   //dispersion parameters for generalised poisson distribution for clutch size
//   //PARAMETER(poisson_disper); 
//   
//   int nObs  = LayDateApril.size();
//   int nYear = NLEVELS(t);
//   
//   Type nll=Type(0);
//   //through if statement, we turn on or off alpha, theta or omega
//   if (CppAD::Variable(mu_alpha(0)) || CppAD::Variable(mu_theta(0)) || CppAD::Variable(mu_omega(0)))
//   {
//     int d = 3;
//     
//     UNSTRUCTURED_CORR_t<Type> neg_log_density(parCor);
//     // Sigma is the correlation matrix of the white noises for VAR1 mu_alpha, mu_theta and mu_omega
//     matrix<Type> Sigma=neg_log_density.cov(); 
//     
//     vector<Type> vecSigma(d*d);
//     
//     for (int i=0; i<d; i++)
//       for (int j=0; j<d; j++) {
//         vecSigma(i+j*d) = Sigma(i,j);}
//       
//       ADREPORT(vecSigma);
//     ADREPORT(Sigma);
//     
//     matrix<Type> A(d*d,d*d);
//     for (int i=0; i<d; i++)
//       for (int j=0; j<d; j++)
//         for (int k=0; k<d; k++)
//           for (int l=0; l<d; l++)
//             A(i*d+k, j*d+l) = -Phi(i,j)*Phi(k,l);
//     for (int i=0; i<d*d; i++)
//       A(i,i) += 1;
//     matrix<Type> Ainv = A.inverse();
//     
//     //Gamma0 is the stationary covariance matrix for VAR1 mu_alpha, mu_theta and mu_omega
//     vector<Type> vecGamma0 = Ainv*vecSigma;
//     matrix<Type> Gamma0(d,d);
//     for (int i=0; i<d; i++)
//       for (int j=0; j<d; j++){
//         Gamma0(i,j) = vecGamma0(i+j*d);}
//       ADREPORT(Gamma0);
//     matrix<Type> stdGamma0=sqrt(Gamma0.array()); // report standard deviation instead of variance for mu_alpha, mu_theta and mu_omega
//     ADREPORT(stdGamma0);
//     
//     matrix<Type> mu(nYear,d);
//     for (int t=0; t<nYear; t++){
//       mu(t,0)= mu_alpha(t);
//       mu(t,1)= mu_theta(t);
//       mu(t,2)= mu_omega(t);
//     }
//     
//     // Rescaling such such that mu_alpha, mu_theta and mu_omega are standard normal  
//     vector<Type> sds(3);
//     sds(0) = 1/sqrt(Gamma0(0,0));
//     sds(1) = 1/sqrt(Gamma0(1,1));
//     sds(2) = 1/sqrt(Gamma0(2,2));
//     
//     nll += VECSCALE(MVNORM(Gamma0),sds)(mu.row(0));
//     
//     for (int t=1; t<mu.rows(); t++) {
//       vector<Type> mut = mu.row(t);
//       vector<Type> PhiMu = mu.row(t-1)*Phi;
//       vector<Type> resid = mut - PhiMu;
//       
//       nll += VECSCALE(MVNORM(Sigma),sds)(resid); 
//     }
//   }
//   
//   //random year effect for zip contribution to likelihood
//   if (CppAD::Variable(mu_zip(0)))
//     nll-= sum(dnorm(mu_zip,Type(0),Type(1),1)); //iid N(0,1) random year effect
//   
//   
//   //random mother effect contribution to likelihood
//   //iid N(0,1) female id effects
//   if (CppAD::Variable(epsi_surv1(0)))
//     nll-= sum(dnorm(epsi_surv1,Type(0),Type(1),1)); 
//   if (CppAD::Variable(epsi_surv2(0)))
//     nll-= sum(dnorm(epsi_surv2,Type(0),Type(1),1)); 
//   if (CppAD::Variable(epsi_surv3(0)))
//     nll-= sum(dnorm(epsi_surv3,Type(0),Type(1),1)); 
//   if (CppAD::Variable(epsi_mothersurv(0)))
//     nll-= sum(dnorm(epsi_mothersurv,Type(0),Type(1),1)); 
//   
//   //contribution of random intercepts and slopes of directional selection to likelihood
//   // vector<Type> alphaRan(2);
//   // matrix<Type> CorRan(2,2);  
//   // CorRan(0,0) = 1; CorRan(0,1) = rhoRan;
//   // CorRan(1,0) = rhoRan; CorRan(1,1) = 1;
//   // for (int t=0; t<nYear; ++t){
//   //   if (CppAD::Variable(logSdRan(0))){
//   //     alphaRan(0)=alpha0_t(t);
//   //     alphaRan(1)=alpha1_t(t);
//   //     nll += VECSCALE(MVNORM(CorRan),sdVecRan)(alphaRan); }
//   // }
//   
//   
//   matrix<Type> survival(nObs,4); 
//   matrix<Type> zip(nObs,3);   
//   
//   for (int s=0; s<4; s++){
//     for (int i=0; i<nObs; i++){
//       
//       Type etaAlpha = alpha_0(s)      
//       +alpha_BCI2(s)*BCI2(t(i))    
//       +alpha_BCI3(s)*BCI3(t(i))      
//       +alpha_dens(s)*dens(t(i))
//       +alpha_clusize(s)*y(i,0);
//       
//       if (CppAD::Variable(logalpha_mother(s))){
//         if (s==0) {etaAlpha += alpha_mother(s)*epsi_surv1(j(i));} 
//         else if (s==1) {etaAlpha += alpha_mother(s)*epsi_surv2(j(i));}
//         else if (s==2) {etaAlpha += alpha_mother(s)*epsi_surv3(j(i));}
//         else {etaAlpha += alpha_mother(s)*epsi_mothersurv(j(i));} 
//       }
//       
//       
//       if (CppAD::Variable(logalpha_year(s)))
//         etaAlpha += alpha_year(s)*mu_alpha(t(i)); 
//       
//       Type eta=etaAlpha;
//       
//       Type etaTheta = theta_0(s)+theta_temp(s)*temp(t(i))+theta_dens(s)*dens(t(i));
//       
//       if (CppAD::Variable(theta_peakd(s)))
//         etaTheta +=  theta_peakd(s)*peakd(t(i));
//       if (CppAD::Variable(theta_clusize(s)))
//         etaTheta += theta_clusize(s)*y(i,0);  
//       if (CppAD::Variable(logtheta_year(s)))
//         etaTheta += theta_year(s)*mu_theta(t(i)); 
//       
//       Type etaOmega = omega_0(s)+omega_temp(s)*temp(t(i));
//       
//       if (CppAD::Variable(logomega_year(s)))
//         etaOmega += omega_year(s)*mu_omega(t(i)); 
//       
//       
//       // if (s==0) { //for poission-distributed clutch size
//       //   if (CppAD::Variable(theta_0(s))){
//       //     eta -= pow((LayDateApril(i)-etaTheta)/exp(etaOmega),2)/2;}else if
//       //       (CppAD::Variable(alpha1_laydate)){
//       //       eta += alpha1_laydate*LayDateApril(i);
//       //       if(CppAD::Variable(alpha0_t(0)))  eta += alpha0_t(t(i));
//       //       if(CppAD::Variable(alpha1_t(0)))  eta += alpha1_t(t(i))*LayDateApril(i);
//       //     }
//       // 
//       //     
//       //     poi_lamb(i) = exp(eta);
//       //     nll-= dgpois(y(i,0), poi_lamb(i),poisson_disper, true);}else{ 
//       
//       // for survival probabilitis including mother survival
//       eta += offset(s); 
//       if (CppAD::Variable(theta_0(s))){
//         eta += pow((LayDateApril(i)-etaTheta)/exp(etaOmega),2)/2; }
//       
//       
//       //cumulative hazard function
//       Type cuHazard = exp(eta); 
//       survival(i,s) = exp(-cuHazard);
//       
//       // 0 observations are removed since they do not cotribute to the probability likelihood
//       //alpha and beta are shape parameters in beta function
//       //disPar are dispersion parameters 
//       if (s<3 && y(i,s) > 0){
//         Type alpha = survival(i,s)*(y(i,s) - 1)/(disPar(s) - 1) + 1e-100;
//         Type beta  = (1 - survival(i,s))*(y(i,s) - 1)/(disPar(s) - 1) + 1e-100;
//         
//         if (CppAD::Variable(zip_0(s))){
//           // logit of zero-inflated probability
//           Type  zip_fun = zip_0(s)           
//           +beta_RLD(s)*RLD(i)      
//           +beta_dens(s)*dens(t(i)) 
//           +beta_diffSurvFem(s)*diffSurvFem(t(i))
//           +beta_peakd(s)*peakd(t(i))     
//           +beta_LayDate(s)*LayDateApril(i)
//           +offset(s);   
//           
//           //clutch size effect on zip
//           if (CppAD::Variable(beta_clusize(s)))
//             zip_fun += beta_clusize(s)*y(i,0);  
//           
//           
//           //random year effect on zip
//           if (CppAD::Variable(logbeta_year(s)))
//             zip_fun += beta_year(s)*mu_zip(t(i)); 
//           
//           zip(i,s) = 1-exp(-exp(zip_fun)); 
//           
//           
//           
//           nll -= dzibetabinom(y(i,s+1), y(i,s), alpha, beta, zip(i,s),true);}
//         else{ nll -= dbetabinom(y(i,s+1), y(i,s), alpha, beta,true);}
//       } 
//       
//       else if (s==3 && brood(i)==0){
//         nll -= dbinom(y(i,4),Type(1), survival(i,3), true);}
//       
//       
//     }
//   }
//   
//   // ADREPORT(poi_lamb);
//   //  ADREPORT(survival);
//   // ADREPORT(zip);
//   return nll;
// }
// 
//                           
//                           
//                           
//                           