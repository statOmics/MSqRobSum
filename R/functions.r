##' Robust differential abundance analysis for label-free quantitative proteomics and robust peptide abundance summarization
##' 
##' Perform robust differential abundance analysis directly on peptide abundances from label-free quantitative proteomics experiment. Or first do a robust summarization on these peptide abundances to protein abundances and perform a robust differential abundance analysis on these summarized protein abundances.
##'
#' @param data MSnset object or dataframe or with at least folowing 3 columns:
#'\describe{
#' \item{expression}{Expression values in log scale.}
#' \item{sample}{Which sample/run the measurement belongs to. Only required when mode = ``sum'' or ``msqrobsum''.}
#' \item{feature}{Which feature (eg. peptide id) the measurement belongs to. Only required when mode = ``sum'' or ``msqrobsum''.}
#' \item{Any variables specified in formulas.}{}
#' }
#' @param formulas Vector of formulas.
#' These are the msqrob model specifications.
#' A two-sided linear ``lme4'' formula object describing both the fixed-effects and random-effects part of the model, with the response on the left of a ~ operator and the terms, separated by + operators, on the right.
#' Random-effects terms are distinguished by vertical bars (|) separating expressions for design matrices from grouping factors (eg. (1 | treatment)). See ``lme4'' package for more details.
#' When multiple models are specified then the first formula is tried first to fit the data.
#' When this fails, the second is tried, etc.
#' @param group_vars Character vector of variable names.
#' The variables used to group the data (eg. protein id). A model will be fitted for each group.
#' @param contrasts Numeric matrix with contrasts or character with variable name.
#' When a variable name is specified then this should correspond to categorial variable (eg. treatment) an should be specified in the model as a random effect.
#' Every possible contrast will then be calculated.
#' The contrast matrix should also only involve categorial parameters specified as random effect.
#' This is because the reference level in the model can change between groups (eg. proteins) due to missing category levels
#' @param mode Character.
#' ``'msqrobsum''' Summarization and MSqRob analysis is performed on the data.
#' ``'msqrob''' Only MSqRob analysis is performed on the data.
#' ``'sum''' Only Summarization is performed on the data.
#' @param robust_lmer_iter Integer or ``'auto'''. Number of iterations used for robust estimation in MSqRob (M-estimation with Huber weights).
#' when set to ``'auto''', defaults to 1 if ``mode = msqrobsum'' and 20 if ``mode = msqrob''
#' @param squeeze_variance Logical. ``TRUE'' if you want to squeeze the residual standard deviation of all models should be squeezed towards a common value
#' @param p_adjust_method Character. Correction method for multiple testing.
#' Defaults to "fdr". See ``fdrtool::p.adjust'' for more information an all available methods.
#' @param keep_model Logical. ``TRUE'' (default) if you want to keep all lme4 models in the output. (memory heavy)
#' @param rlm_args Named list. All parameters to be passed to the `rlm` function used in the summarization step. Default parameters when empty list. See ``MASS::rlm'' for more information on all parameters and default settings.
#' @param lmer_args Named list. All parameters to be passed to the `lmer` function used in the MSqRob analysis. Default parameters when empty list. See ``lme4::lmer'' for more information on all parameters and default settings.
#' @param parallel_args Named list. All parameters to be passed to the `plan` function from the `future` package which allows for parallelization. Set ``strategy = 'multisession'' to allow parallelization using all available cores (default). Set the ``workers'' parameter to an integer to choose the number of cores to be used. Set ``strategy = sequential'' to disable parallelization. See ``future::plan'' for more information on all available perallelization strategies and other parameters with their default settings.
#' @return A data frame.
#' Following columns are present:
#' \describe{
#'\item{group vars}{}
#'\item{data}{}
#'\item{data_summarized}{}
#'\item{formula}{}
#'\item{df}{}
#'\item{sigma}{}
#'\item{intercept}{}
#'\item{sigma_post}{}
#'\item{df_prior}{}
#'\item{sigma_prior}{}
#'\item{df_post}{}
#'\item{contrasts}{}
#' }
#' @export
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @import lme4
#' @import future.apply
#' @import future
#' @import limma
#'
#' @examples
#' ## Robust summarization from peptide intensities to protein summaries on
##' ## the build-in data set with peptide intensities from 100 proteins.
##' ## For only 100 proteins we will not benefit from parallezation because
##' ## robust summarization is a fairly fast routine.
##' results1 <- msqrobsum( data = peptide_intensities
##'                     , mode = 'sum'
##'                     , group_vars = 'protein'
##'                     , parallel_args = list(strategy = 'sequential'))
##' 
##' ## MSqRobSum analysis
##' ## There are 20 samples belonging to 5 different conditions.
##' ## Differential expression is tested between all conditions.
##' form = expression ~ (1|condition)
##' results2 <- msqrobsum(data = peptide_intensities
##'                     , formulas = form
##'                     , mode = 'msqrobsum'
##'                     , group_vars = 'protein'
##'                     , contrasts = 'condition'
##'                     , parallel_args = list(strategy = 'sequential'))
##' 
##' ## MSqRob analysis
##' ## There are 20 samples belonging to 5 different conditions.
##' ## Since there is no prior summarization from peptide to protein intensities.
##' ## The model has to take into account the sample and feature (peptide) effects
##' form =  c(expression ~ (1|condition) + (1|sample) + (1|feature), expression ~ (1|condition))
##' ## Differential expression is tested between all conditions.
##' ## Fitting the full MSqRob models takes longer then the simplified models in MSqRobSum.
##' ## Therefore it's suggested that you allow for parallelization,
##' ## especially if you have big data sets with many samples and thousands of proteins.
##' ## eg. if you have 2 available processing cores.
##' results3 <- msqrobsum(data = peptide_intensities
##'                     , formulas = form
##'                     , mode = 'msqrob'
##'                     , group_vars = 'protein'
##'                     , contrasts = 'condition'
##'                     , parallel_args = list(strategy = 'multisession', workers = 2))
msqrobsum <- function(
                      data, formulas, group_vars = 'protein', contrasts = NULL
                    , mode = c('msqrobsum','msqrob', 'sum')
                    , robust_lmer_iter = 'auto'
                    , squeeze_variance = TRUE
                    , p_adjust_method = c("BH",p.adjust.methods)
                    , keep_model = FALSE
                    , rlm_args = list(maxit = 20L)
                    , lmer_args = list(control = lmerControl(calc.derivs = FALSE))
                    , parallel_args = list(strategy = 'multisession')
                      ## EXPERIMENTAL OPTIONS
                    , type_df = 'traceHat'
                      ## , type_df = c('traceHat','conservative')
                    , squeeze_covariate = FALSE
                    , fit_fun = do_mm
                      ){
  rlang::exec(plan, !!!parallel_args)
  mode = match.arg(mode)

  if (robust_lmer_iter == 'auto')
    robust_lmer_iter <- ifelse(mode == 'msqrob', 20L, 1L)
  type_df = match.arg(type_df)
  p_adjust_method = match.arg(p_adjust_method)
  if (missing(formulas) & mode == 'sum') formulas = ''
  if ((mode != 'sum') & length(unlist(map(formulas,findbars))) == 0)
    stop('msqrobsum only works with formulas with at least 1 random effect specified')
  if (is(data,'MSnSet')) data <- MSnSet2df(data)
  if (is(contrasts, 'character')){
    if (!all(c(formulas) %>% map_lgl(~{contrasts %in% map_chr(findbars(.),all.vars)})))
      stop('Contrasts can only be constructed from variable name when the variable is specified as random effect')
    contrasts <- make_simple_contrast(data, contrasts)}
  group_vars <- map_chr(group_vars, ~match.arg(.,names(data)))
  model_vars <- c(formulas) %>% map(all.vars) %>% unlist %>% unique

  ## select only necessary columns needed to save memory
  df <- data %>% select(group_vars,model_vars, 'expression', 'feature','sample') %>% 
    group_by_at(group_vars) %>% nest %>%
    mutate(mm = future_lapply(data, fit_fun, robust_lmer_iter = robust_lmer_iter
                            , rlm_args = rlm_args, lmer_args =lmer_args
                            , contrasts = contrasts, keep_model = keep_model
                            , formulas = formulas, mode = mode)) %>%
    unnest(mm) %>% ungroup
  if (mode == 'sum') return(df)
  ## Return also failed ones afterward
  df_failed <- filter(df, is.na(df))
  df <- filter(df, !is.na(df))
  if(!nrow(df)) {warning("No models could be fitted"); return(df_prot_failed)}
  ## Squeeze variance
  df <- mutate(df, sigma_post = sigma, df_prior = 0, sigma_prior = 0)
  if(squeeze_variance) {
    ## no shrinkage if df < 1 because greatly influences emp. bay.
    ## TODO check only works when at least 3 protein
    id = df$df >= 1L
    sq = squeezeVar(df$sigma[id]^2, df$df[id])
    ## experimental: use intercept of model as covariate to
    ## squeeze variance to correct for  mean variance
    if(squeeze_covariate) sq = squeezeVar(df$sigma[id]^2, df$df[id], covariate = df$intercept[id])
    df[id,] = mutate(df[id,]
                   , sigma_prior = sqrt(sq$var.prior), sigma_post = sqrt(sq$var.post), df_prior = sq$df.prior)}
  if(type_df == "conservative"){
    ## Calculate df on protein level, assumption is that there is only one protein value/run,
    ## TODO remove conservative option? Now broken, fix vars = colnames(...)
    df <- mutate(df, df_prior = 0
               , df = map2_dbl(data, model,~calculate_df(.x,.y, vars = colnames(pData(msnset)))))}
  df <- mutate(df, df_post = df + df_prior)
  
  ## Calculate q and p values for contrasts (if contrasts are found)
  df <- df %>% select(group_vars ,df_post,sigma_post, contrasts) %>% unnest %>% 
    ## skip if there are no contrasts
    {if(nrow(.)) {
       mutate(., se = sigma_contrast * sigma_post
            , t = logFC/se
            , pvalue = pt(-abs(t), df_post) * 2) %>%
         select(-sigma_post, - df_post) %>%
         group_by(contrast) %>%
         mutate(qvalue = p.adjust(pvalue, method = p_adjust_method)) %>%
         group_by_at(group_vars) %>% nest(.key = contrasts)
     } else .} %>%
    select(group_vars, contrasts) %>%
    left_join(select(df, -contrasts),., by = group_vars)
  
  ## mutate(df_failed, contrasts = list(tibble())) %>% bind_rows(df,.)
  ## give empty tibble instead of NULL when no contrast (nicer to work with posthoc)
  bind_rows(df,df_failed) %>%
    mutate(contrasts = map(contrasts, ~{if (is(.x,'tbl')) .x else tibble()}))
}

#' @export
MSnSet2df <- function(msnset){
  ## Converts Msnset to a tidy dataframe
  ## Always creates feature and vector column so these shouldn't be defined by user.
  ## convenient for downstream analysis steps.
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Package \"Biobase\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if(any(c("sample", "feature", "expression") %in% c(colnames(fData(msnset)),colnames(pData(msnset))))){
    stop("Column names in the \"fData\" or \"pData\" slot of the \"msnset\" object cannot be named
         \"sample\", \"feature\" or \"expression\". Please rename.")
  }

  dt <- as.data.frame(Biobase::exprs(msnset)) %>% mutate(feature = rownames(.)) %>%
    gather(sample, expression, - feature, na.rm=TRUE)
  dt <- fData(msnset) %>% mutate(feature = (rownames(.))) %>% left_join(dt,. , by = 'feature')
  dt <- pData(msnset) %>% mutate(sample = rownames(.)) %>% left_join(dt,. , by = 'sample')
  as_data_frame(dt)
}

#' @import lme4
#' @import tibble
do_mm <- function(d, formulas, contrasts, mode ='msqrobsum', robust_lmer_iter = 1L
                , rlm_args = list(), lmer_args = list(), keep_model = TRUE){
  out <- list()
  if (mode != 'msqrob') {
    d <- mutate(d,expression = rlang::exec(summarizeRobust, expression, feature, sample, !!!rlm_args)) %>%
      select(-feature) %>% distinct
    out$data_summarized <- list(d)
  }
  if (mode == 'sum') return(new_tibble(out ,nrow = 1L))
  for (form in c(formulas)){
    ## TODO do lm fit if there are no random effects (all output inside loop, do while loop)
    model <- try(do_lmerfit(d, form,  robust_lmer_iter,lmer_args), silent = TRUE)
    if (is(model,"lmerMod")) break #else message(paste0('Cannot fit ', format(form)))
  }
  if (keep_model) out$model <- list(model)
  if (!is(model,"lmerMod")) return(new_tibble(out ,nrow = 1L))
  out$formula <- as.character(list(attributes(model@frame)$formula))
  out$df <- getDf(model)
  out$sigma <- sigma(model)
  out$contrasts <- list(get_contrasts(model,contrasts))
  out$intercept <- model@beta[1]
  new_tibble(out , nrow = 1L) 
}

#' @importFrom MASS psi.huber
#' @import lme4
do_lmerfit <- function(df, form, robust_lmer_iter, args_lmer = list()){
  tol <- 1e-6
  fit <- rlang::exec(lmer,form, data = df, !!!args_lmer)
  ##Initialize SSE
  sseOld <- fit@devcomp$cmp['pwrss']
  while (robust_lmer_iter > 0){
    robust_lmer_iter= robust_lmer_iter-1
    res <- resid(fit)
    fit@frame$`(weights)` <- psi.huber(res/(mad(res,0)))
    fit <- refit(fit)
    sse <- fit@devcomp$cmp['pwrss']
    if(abs(sseOld-sse)/sseOld <= tol) break
    sseOld <- sse
  }
  fit
}

#' @importFrom MASS rlm
summarizeRobust <- function(expression, feature, sample, ...) {
  ## Assumes that intensities mx are already log-transformed
  ## characters are faster to construct and work with then factors
  feature <- as.character(feature)
  ##If there is only one 1 peptide for all samples return expression of that peptide
  if (length(unique(feature)) == 1L) return(expression)
  sample <- as.character(sample)
  ## modelmatrix breaks on factors with 1 level so make vector of ones (will be intercept)
  if (length(unique(sample)) == 1L) sample <- rep(1,length(sample))

  ## sum contrast on peptide level so sample effect will be mean over all peptides instead of reference level
  X <- model.matrix(~ -1 + sample + feature,contrasts.arg = list(feature = 'contr.sum'))
  ## MASS::rlm breaks on singular values.
  ## check with base lm if singular values are present.
  ## if so, these coefficients will be zero, remove this column from model matrix
  ## rinse and repeat on reduced modelmatrix untill no singular values are present
  repeat {
    fit <- .lm.fit(X,expression)
    id <- fit$coefficients != 0
    X <- X[ , id, drop =FALSE]
    if (!any(!id)) break
  }
  ## Last step is always rlm
  fit <- rlm(X, expression, ...)
  ## Only return the estimated effects effects as summarised values
  ## sampleid = seq_along(unique(sample))
  ## return(X[,sampleid,drop = FALSE] %*% fit$coefficients[sampleid])
  fit$coefficients[paste0('sample',sample)]
}

#' @import tibble
get_contrasts <- function(model, contrasts){
  ## TODO allow for contrasts from fixed effects
  ## tricky because reference level can change
  betaB <- getBetaB(model)
  vcov <- getVcovBetaBUnscaled(model)
  coefficients <- names(betaB)
  id <- coefficients %in% rownames(contrasts)
  if (!any(id)) return(new_tibble(list(), nrow = 0))
  coefficients <- coefficients[id]
  vcov <- vcov[id,id]
  betaB <- betaB[id]

  ## check for which contrasts I have data
  missing_levels <- !(rownames(contrasts) %in% coefficients)
  id <- !apply(contrasts,2,function(x){any(x[missing_levels] != 0)})
  contrasts <- contrasts[coefficients, id, drop = FALSE]
  ## If no contrasts could be found, terminate
  if (is.null(colnames(contrasts))) return(new_tibble(list(), nrow = 0))
  new_tibble(list(contrast = colnames(contrasts)
                , logFC = logFC <- (t(contrasts)%*%betaB)[,1]
                , sigma_contrast = sqrt(diag(t(contrasts)%*%vcov%*%contrasts))
                  ),nrow = sum(id))
}

#' @import purrr
getBetaB <- function(model) {
  betaB <- c(as.vector(getME(model,"beta")),as.vector(getME(model,"b")))
  names(betaB) <- c(colnames(model@pp$X), unlist(imap(model@flist,~{paste0(.y,levels(.x))})))
  betaB
}

#' @import lme4
getVcovBetaBUnscaled <- function(model){
  X <- getME(model,"X")
  Z <- getME(model,"Z")
  vcovInv <- Matrix::crossprod(cbind2(X,Z))
  Ginv <- Matrix::solve(Matrix::crossprod(getME(model,"Lambda")) +
                        Matrix::Diagonal(ncol(Z),1e-18))
  i <- -seq_len(ncol(X))
  vcovInv[i,i] <- vcovInv[i,i]+Ginv
  as.matrix(Matrix::solve(vcovInv))
}

calculate_df <- function(df, model, vars){
  ## Get all the variables in the formula that are not defined in vars
  form <- attributes(model@frame)$formula
  vars_formula <- all.vars(form)
  vars_drop <- vars_formula[!vars_formula %in% vars]
  ## Sum of number of columns -1 of Zt mtrix of each random effect that does not involve a variable in vars_drop
  mq <- getME(model,'q_i')
  id <- !map_lgl(names(mq),~{any(stringr::str_detect(.x,vars_drop))})
  p <- sum(mq[id]) - sum(id)
  ## Sum of fixed effect parameters that do not involve a variable in vars_drop
  mx <- getME(model,'X')
  id <- !map_lgl(colnames(mx),~{any(stringr::str_detect(.x,vars_drop))})
  p <- p + sum(id)

  ## n is number of sample because 1 protein defined per sample
  n <- n_distinct(df$sample)
  n-p
}

#' @import lme4
getDf <- function(object){
  w <- object@frame$"(weights)"
  if (is.null(w)) w <- 1
  sigma <- sigma(object)
  sum((resid(object)* sqrt(w))^2)/sigma^2
}

#' @export
make_simple_contrast <- function(data, contrast_var){
  c <- pull(data,contrast_var) %>% unique %>% paste0(contrast_var, .) %>% sort %>% as.factor
  comp <- combn(c,2,simplify = FALSE)
  condIds <- map(comp, ~which(c %in% .x))
  L <- rep(0,nlevels(c))
  L <- sapply(comp,function(x){L[x]=c(-1,1);L})
  rownames(L) <- levels(c)
  colnames(L) <- map_chr(comp, ~paste(rev(.x),collapse = '-'))
  L
}

##################
## Experimental ##
##################
## fit with inverse variances of sample as weights

#' @import lme4
do_mm2 = function(d, formulas, contrasts, mode ='msqrobsum', robust_lmer_iter = 1L
                , rlm_args = list(), lmer_args = list(), keep_model = TRUE){
  out = list()
  ## print('\\n1\\n')
  if (mode != 'msqrob') {
    d = select(d,-feature,-expression) %>%
      bind_cols(summarizeRobust2(d$expression, d$feature, d$sample)) %>% 
      group_by_at(vars(-expression,-residuals)) %>%
      ## summarise(expression = first(expression), weights = 1/(mad(residuals,0)^2)) %>%
      summarise(expression = first(expression), weights = 1/(mean(residuals^2))) %>%
      ungroup #%>%#glimpse %>%
    ## mutate(weights = ifelse(is.infinite(weights), min(weights),weights)
    ##      , weights = ifelse(weights > 10e9 , min(weights), weights)
    ##      , weights = weights/max(weights) # sometimes all residuals are very low
    ##      , weights = ifelse(is.na(weights), 1, weights) # models have a NA coef, no residuals are calculated
    ## ) # %>% glimpse
    if (any(is.infinite(d$weights) |is.na(d$weights) | d$weights > 10e6)) d$weights = 1
    ## d$weights = d$weights/max(d$weights)
    out$data_summarized = list(d)
  }
  if (mode == 'sum') return(new_tibble(out ,nrow = 1L))
  ## print('\\n2\\n')
  for (form in c(formulas)){
    model = try(do_lmerfit2(d, form,  robust_lmer_iter,lmer_args), silent = TRUE)
    if (is(model,"lmerMod"))  break else message(paste0('Cannot fit ', format(form)))
  }
  if (keep_model) out$model = list(model)
  if (!is(model,"lmerMod")) return(new_tibble(out ,nrow = 1L))
  out$formula = as.character(list(attributes(model@frame)$formula))
  out$df = getDf(model)
  out$sigma = sigma(model)
  out$contrasts = list(get_contrasts(model,contrasts))
  out$intercept = model@beta[1]
  new_tibble(out , nrow = 1L) 
}

#' @importFrom MASS psi.huber
#' @import lme4
do_lmerfit2 = function(df, form, robust_lmer_iter, args_lmer = list()){
  tol = 1e-6
  ## print(glimpse(df))
  ## fit <- rlang::exec(lmer,form, data = df, !!!args_lmer)
  ## print('\\n3\\n')
  fit <- lmer(formula = form, data = df, weights = weights)
  ## print('\\n4\\n')
  ##Initialize SSE
  res <- resid(fit)
  sseOld <- fit@devcomp$cmp['pwrss']
  while (robust_lmer_iter > 0){
    robust_lmer_iter= robust_lmer_iter-1
    fit@frame$`(weights)` <- df$weights * psi.huber(res/(mad(res)))
    fit <- refit(fit)
    res <- resid(fit)
    sse <- fit@devcomp$cmp['pwrss']
    if(abs(sseOld-sse)/sseOld <= tol) break
    sseOld <- sse
  }
  fit
}

#' @importFrom MASS rlm
summarizeRobust2 = function(expression, feature, sample, ...) {
  ## Assumes that intensities mx are already log-transformed
  ## characters are faster to construct and work with then factors
  feature <- as.character(feature)
  ##If there is only one 1 peptide for all samples return expression of that peptide
  if (length(unique(feature)) == 1L) return(list(expression = expression, residuals = rep(1, length(expression))))
  sample <- as.character(sample)
  ## modelmatrix breaks on factors with 1 level so make vector of ones (will be intercept)
  if (length(unique(sample)) == 1L) sample = rep(1,length(sample))

  ## sum contrast on peptide level so sample effect will be mean over all peptides instead of reference level
  X = model.matrix(~ -1 + sample + feature,contrasts.arg = list(feature = 'contr.sum'))
  ## MASS::rlm breaks on singular values.
  ## check with base lm if singular values are present.
  ## if so, these coefficients will be zero, remove this column from model matrix
  ## rinse and repeat on reduced modelmatrix untill no singular values are present
  repeat {
    fit = .lm.fit(X,expression)
    id = fit$coefficients != 0
    X = X[ , id, drop =FALSE]
    if (!any(!id)) break
  }
  ## Last step is always rlm
  fit = rlm(X, expression, ...)
  ## Only return the estimated effects effects as summarised values
  ## sampleid = seq_along(unique(sample))
  ## return(X[,sampleid,drop = FALSE] %*% fit$coefficients[sampleid])
  ## print(fit)
  ## print(resid(fit))
  list(expression = as.vector(fit$coefficients[paste0('sample',sample)]),
       residuals = fit$residuals
       )
}

#################################################
#' @import lme4
do_mm3 = function(d, formulas, contrasts, mode ='msqrobsum', robust_lmer_iter = 1L
                , rlm_args = list(), lmer_args = list(), keep_model = TRUE){
  out = list()
  if (mode != 'msqrob') {
    d = mutate(d,expression = rlang::exec(summarizeRobust, expression, feature, sample,!!!rlm_args)) %>%
      select(-feature) %>% distinct
    out$data_summarized = list(d)
  }
  if (mode == 'sum') return(new_tibble(out ,nrow = 1L))
  for (form in c(formulas)){
    model = try(do_lmerfit3(d, form,  robust_lmer_iter,lmer_args), silent = TRUE)
    if (is(model,"lmerMod"))  break else message(paste0('Cannot fit ', format(form)))
  }
  if (keep_model) out$model = list(model)
  if (!is(model,"lmerMod")) return(new_tibble(out ,nrow = 1L))
  out$formula = as.character(list(attributes(model@frame)$formula))
  out$df = getDf(model)
  out$sigma = sigma(model)
  out$contrasts = list(get_contrasts(model,contrasts))
  out$intercept = model@beta[1]
  new_tibble(out , nrow = 1L) 
}

#' @importFrom MASS psi.huber
#' @import lme4
do_lmerfit3 = function(df, form, robust_lmer_iter, args_lmer = list()){
  tol = 1e-6
  fit <- rlang::exec(lmer,form, data = df, !!!args_lmer)
  ##Initialize SSE
  sseOld <- fit@devcomp$cmp['pwrss']
  while (robust_lmer_iter > 0){
    robust_lmer_iter= robust_lmer_iter-1
    res <- resid(fit)
    fit@frame$`(weights)` <- psi.huber(res/(mad(res,0)))
    fit <- refit(fit)
    sse <- fit@devcomp$cmp['pwrss']
    if(abs(sseOld-sse)/sseOld <= tol) break
    sseOld <- sse
  }
  fit
}
