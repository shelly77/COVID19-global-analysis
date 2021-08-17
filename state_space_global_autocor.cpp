
// This c++ file is self-defined model in TMB format


#include <TMB.hpp>

// we comment out these self-defined distribution functions, they are not used in this project
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
//template<class Type> // define beta-binomial distribution
//inline Type dbetabinom(const Type &x, const Type &size, const Type &alpha, const Type &beta,
//                       int give_log=0) // size is the total number of trials and x is the number of success
//{
// Type logres =
//   lfactorial(size)
// - lfactorial(x)
// - lfactorial(size - x)
//  + lgamma(alpha + x)
//  + lgamma(alpha + beta)
//  - lgamma(alpha + beta + size)
//  - lgamma(alpha)
//  - lgamma(beta)
//  + lgamma(beta + size - x+ 1e-100);
//  if (give_log) return logres; else return exp(logres);
//}
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
  
  // data in TMB format
  DATA_VECTOR(lognewcasessmooth); 
  
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
  DATA_VECTOR(pweekcases);
  DATA_VECTOR(population);
  DATA_VECTOR(dens);
  DATA_VECTOR(mage);
  DATA_VECTOR(gdp); 
  DATA_VECTOR(newtests);
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
  
  // fixed effect parameters
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
  
  PARAMETER(beta_pweekcases);

  PARAMETER(beta_population);
  PARAMETER(beta_dens);
  PARAMETER(beta_mage);

  PARAMETER(beta_gdp);

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
  
  //random country effects (random intercepts)
  PARAMETER_VECTOR(mu_spatioLevel);
  PARAMETER(logsd_spatioLevel);   Type sd_spatioLevel = exp(logsd_spatioLevel);  ADREPORT(exp(2*logsd_spatioLevel));
  
  //random slopes of days
  PARAMETER_VECTOR(mu_spatioLevel_days);
  PARAMETER(logsd_randm_slop_days);   Type sd_randm_slop_days = exp(logsd_randm_slop_days);  ADREPORT(exp(2*logsd_randm_slop_days));
  
  //random slopes of squared days
  PARAMETER_VECTOR(mu_spatioLevel_sqdays);
  PARAMETER(logsd_randm_slop_sqdays);   Type sd_randm_slop_sqdays = exp(logsd_randm_slop_sqdays); ADREPORT(exp(2*logsd_randm_slop_sqdays));
  
  //random slopes on mobility
  PARAMETER_VECTOR(mu_spatioLevel_mobility); ADREPORT(beta_mobility+mu_spatioLevel_mobility);
  PARAMETER(logsd_randm_slop_mobility);   Type sd_randm_slop_mobility= exp(logsd_randm_slop_mobility); ADREPORT(exp(2*logsd_randm_slop_mobility));
  
  //random slopes of one week lagged mobility
  PARAMETER_VECTOR(mu_spatioLevel_mobility_lag1); ADREPORT(beta_mobility_lag1+mu_spatioLevel_mobility_lag1);
  PARAMETER(logsd_randm_slop_mobility_lag1);   Type sd_randm_slop_mobility_lag1 = exp(logsd_randm_slop_mobility_lag1); ADREPORT(exp(2*logsd_randm_slop_mobility_lag1));
  
  //random slopes of two weeks lagged mobility
  PARAMETER_VECTOR(mu_spatioLevel_mobility_lag2); ADREPORT(beta_mobility_lag2+mu_spatioLevel_mobility_lag2);
  PARAMETER(logsd_randm_slop_mobility_lag2);   Type sd_randm_slop_mobility_lag2 = exp(logsd_randm_slop_mobility_lag2); ADREPORT(exp(2*logsd_randm_slop_mobility_lag2));
  
  //random slopes on temperature
  PARAMETER_VECTOR(mu_spatioLevel_temp);
  PARAMETER(logsd_randm_slop_temp);   Type sd_randm_slop_temp = exp(logsd_randm_slop_temp); ADREPORT(exp(2*logsd_randm_slop_temp));
  
  //random slopes on squared temperature
  PARAMETER_VECTOR(mu_spatioLevel_sqtemp);
  PARAMETER(logsd_randm_slop_sqtemp);   Type sd_randm_slop_sqtemp = exp(logsd_randm_slop_sqtemp); ADREPORT(exp(2*logsd_randm_slop_sqtemp));
  
  // parameters in ARCH(2)
  PARAMETER(logalpha0); Type alpha0 = exp(logalpha0); ADREPORT(exp(2*logalpha0));
  PARAMETER(logalpha1); Type alpha1 = exp(logalpha1); ADREPORT(exp(2*logalpha1));
  PARAMETER(logalpha2); Type alpha2 = exp(logalpha2); ADREPORT(exp(2*logalpha2));
  
  int nObs  = lognewcasessmooth.size();
  // int nspatioLevel = NLEVELS(spatioLevel);
  
  Type nll=Type(0);
 
  //contribution of random country effects to likelihood
  if (CppAD::Variable(mu_spatioLevel(0)))
    nll-= sum(dnorm(mu_spatioLevel,Type(0),sd_spatioLevel,1)); //iid N(0,1) random year effect
  
  if (CppAD::Variable(mu_spatioLevel_days(0)))
    nll-= sum(dnorm(mu_spatioLevel_days,Type(0),sd_randm_slop_days,1)); //iid N(0,1) random year effect
  
  if (CppAD::Variable(mu_spatioLevel_sqdays(0)))
    nll-= sum(dnorm(mu_spatioLevel_sqdays,Type(0),sd_randm_slop_sqdays,1)); //iid N(0,1) random year effect
  //  }
  
  
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
    
    if (CppAD::Variable(beta_pweekcases))
      eta += beta_pweekcases*pweekcases(i);
    
    if (CppAD::Variable(beta_population))
      eta += beta_population*population(i);
    
    if (CppAD::Variable(beta_dens))
      eta += beta_dens*dens(i);
    
    if (CppAD::Variable(beta_mage))
      eta += beta_mage*mage(i);
    
    if (CppAD::Variable(beta_gdp))
      eta += beta_gdp*gdp(i);
    
    if (CppAD::Variable(beta_newtests))
      eta += beta_newtests*newtests(i);
    
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
      eta += mu_spatioLevel(spatioLevelfactor(i));
    
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
    
    
    state(i) = eta;
    
    error(i) = lognewcasessmooth(i) - state(i);
    
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

