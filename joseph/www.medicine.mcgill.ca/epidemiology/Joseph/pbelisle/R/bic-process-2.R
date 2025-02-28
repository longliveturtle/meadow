### bic.process, version 1.6 (Released April 2016) #########################################################################


bic.process <- function(bic.out,
                        data,
                        n = nrow(data),
                        report.ci = T, 
                        ci.level = 0.95, 
                        report.or = bic.out.family=="binomial",
                        report.rate.ratio = bic.out.family=="poisson",
                        report.hazard.ratio = bic.fct=="surv",
                        report.r2 = bic.out.family=="gaussian",
                        notes = NULL, 
                        formula = character(0),
                        discard.incomplete.models = T,
                        mixed.interactionterms.factor.type = bic.out.factor.type,
                        recode.absent.var.value.as.na=T, n.events=integer(0))
{
  # NOTES: - This function was validated only for bic.glm outputs when glm.family is either gaussian or binomial,
  #          for bic.surv and bicreg outputs
  #
  #
  # Arguments:
  # ----------
  #
  # bic.out           the output to a bic.glm(), a bic.surv() or bicreg() call
  #
  # data              the independent variables data matrix used in bic.glm(), bic.surv() or bicreg() call;
  #                   it may also include outcome variable. (The actual data are not really used, but their names and levels (if any variable is a factor) will be used to relabel regression parameters)
  #
  # n                 sample size; its default value is ok when bic.glm() was used. When bic.surv() is used,
  #                   you must enter a value for n to get a non-NA value in return list (you could, for example,
  #                   use n=nrow(data) if 'data' was the data matrix and did not include any missing value [NA]).
  #
  # notes             any info that you would like to be attached to this fct's list result 
  #                     (e.g. comment(s) on subjects selection)
  #
  # formula           formula used as an argument in bic.glm() or bic.surv() call
  #
  # discard.incomplete.models  whether or not models with interaction terms present in model but one or more corresponding 
  #                            direct effects not in model should be discarded 
  #                            (if true, many parameters will be recalculated: postmean, postsd, probne0, postprob, condpostmean and condpostsd)
  #
  # recode.absent.var.value.as.na  if TRUE, parameter estimate of a variable not in a given model
  #                                will be displayed as NA (rather than 0.000000), which we find
  #                                adds to output readability

  # outputs from home-made function bic.surv.timedep do not need to be bic.process'ed

  if (!is.null(bic.out$model.head)) stop("bic.surv.timedep outputs do not need to be bic.process'ed.")


  extra <- list()

  reference.value <- function(x){x[1]}
  
  first.row <- function(w)
  {
    nrows <- nrow(w)
    out <- w * row(w)
    out[out==0] <- nrows + 1
    out <- apply(out, 2, min)
    out[out>nrows] <- NA
    names(out) <- colnames(w)
    out
  } # end of first.row
  

  # ----- Try to read formula from $call if not given through 'formula=' option ------------------------------

  # See which bic() was called and whether it was called with a formula or not

  bic.fct <- c('glm', 'surv', 'bicreg') # list of BICs accepted by this fct

  call.cmd <- as.character(bic.out$call)
  called.fct <- call.cmd[1]
  called.fct.parts <- unlist(strsplit(called.fct, '.', fixed=T))
  used.formula <- rev(called.fct.parts)[1] == "formula"

  bic.fct <- intersect(bic.fct, called.fct.parts)
  if (length(bic.fct) == 0) stop("Sorry. This function only accepts outputs from either bic.glm, bic.surv or bicreg")
  bic.ref <- paste('=bic', bic.fct, 'ref(', sep='.')

  if (bic.fct == "bicreg")
  {
    intercept.in.model <- T
    order <- rep(1, ncol(data))
    vars.in.model <- colnames(data)
    
    formula.factors <- diag(1, nrow=length(vars.in.model))
    rownames(formula.factors) <- vars.in.model
    colnames(formula.factors) <- vars.in.model
    
    tmp <- deparse(bic.out$call)
    tmp <- gsub(' ','', tmp)
    tmp <- unlist(strsplit(gsub("[()]", ",", tmp), ","))
    tmp.w <- pmatch("drop.factor.levels", tmp)
    if (is.na(tmp.w))
    {
      bic.out.factor.type <- F
    }
    else
    {
      bic.out.factor.type <- !as.logical(unlist(strsplit(tmp[tmp.w], '='))[2])
    }
    
    bic.out.family <- 'gaussian'
  }
  else
  {
    if (used.formula)
    {
      tmp.formula <- call.cmd[2]
      explicit.formula <- length(grep('~', tmp.formula, fixed=T)) > 0

      # if explicitely available in call, use the formula from call rather than that given by user

      if (explicit.formula) formula <- tmp.formula

      if (length(formula) == 0) stop("Please resubmit with the option 'formula=' given the same value it was given at bic.", bic.fct, " call.")
    }
    else if (length(formula) == 0)
    {
      formula <- '.....y.....~.'
    }

    formula.terms <- terms(as.formula(formula), data=data)
    
    intercept.in.model <- as.logical(attr(formula.terms, "intercept"))
    order <- attr(formula.terms, "order")
    
    vars.in.model <- attr(formula.terms, "term.labels") # list all independent variables in model 
                                                        # (including interaction terms)
    which.quoted.varname <- grep("`", vars.in.model)
    vars.in.model <- gsub("`", "", vars.in.model)
                                                        
    formula.factors <- attr(formula.terms, "factors")
    colnames(formula.factors) <- gsub("`", "", colnames(formula.factors))
    rownames(formula.factors) <- gsub("`", "", rownames(formula.factors))
    
    bic.out.factor.type <- bic.out$factor.type
    bic.out.family <- ifelse(bic.fct=='surv', '_survival_', bic.out$family)
  }
  
  # Read OR from $call
  
  tmp <- deparse(bic.out$call)
  tmp <- gsub(' ','', tmp)
  tmp <- unlist(strsplit(gsub("[()]", ",", tmp), ","))
  tmp.w <- pmatch("OR=", tmp)
  if (is.na(tmp.w))
  {
    OR <- 20 # default bicreg, bic.glm and bic.surv OR value
  }
  else
  {
    OR <- as.numeric(unlist(strsplit(tmp[tmp.w], '='))[2])
  }

  # Read censoring variable from formula

  censoring.var <- NULL
  if (bic.fct == "surv")
  {
    censoring.var <- gsub(".*,",  "", formula, perl=T)
    censoring.var <- gsub("\\).*", "", censoring.var, perl=T)
    censoring.var <- gsub(" ", "", censoring.var)
  }


  # Fetch variable names _____________________________________________________________________________________________


  order1.factors.levels <- lapply(unclass(data), levels) # list of order-1 factors levels
  order1.vars.in.model  <- vars.in.model[order==1]       # list of independent variables
  
  namesx.are.short <- bic.out.factor.type & bic.fct != "bicreg"
  
  # change "X.function" to "function", if a variable in the model was called "function"
  # bic.out.namesx <- bic.out$namesx
  bic.out.namesx <- names(bic.out$probne0)
  tmp <- substr(bic.out.namesx, 1, 2)
  tmp.w <- intersect(which(tmp=="X."), which.quoted.varname)
  if (length(tmp.w))
  {
    tmp.names <- bic.out.namesx[tmp.w]
    tmp.names <- substr(tmp.names, 3, nchar(tmp.names))
    tmp.w.dot <- which(substr(tmp.names,nchar(tmp.names),nchar(tmp.names))==".")
    tmp.names[tmp.w.dot] <- substr(tmp.names[tmp.w.dot],1,nchar(tmp.names[tmp.w.dot])-1)
    tmp.ww <- which(!is.na(match(tmp.names,vars.in.model)))
    
    bic.out.namesx[tmp.w[tmp.ww]] <- tmp.names[tmp.ww]
  }
  

  if (namesx.are.short)
  {
    order1.vars.in.model <- intersect(order1.vars.in.model, bic.out.namesx)
  }
  else
  {
    tmp.order1.levels <- unlist(order1.factors.levels)
    tmp.l <- lapply(order1.factors.levels, length)
    tmp.names <- rep(names(tmp.l), tmp.l)
    tmp.longnames <- paste(tmp.names, tmp.order1.levels, sep='')
    tmp.preferred.longnames <- paste(tmp.names, tmp.order1.levels, sep='=')
    tmp.order1var.in.model <- match(tmp.longnames, bic.out.namesx)
    order1.vars.in.model <- c(intersect(order1.vars.in.model, bic.out.namesx), unique(tmp.names[which(!is.na(tmp.order1var.in.model))]))
    tmp.ordered.varnames <- rownames(formula.factors)
    order1.vars.in.model <- intersect(tmp.ordered.varnames, order1.vars.in.model)
  }

  factors <- formula.factors[order1.vars.in.model,,drop=F] # ignore line concerning outcome variable

  # Determine whether each term contains continuous and/or categorical variables

  is.na1 <- function(x){(length(x) == 1) & all(is.na(x))}
  term.is.cont <- matrix(unlist(lapply(order1.factors.levels[order1.vars.in.model],is.null)), nrow=1)

  any.cont  <- as.vector(term.is.cont %*% factors) > 0
  n.categ   <- as.vector((!term.is.cont) %*% factors)
  any.categ <- n.categ > 0

  any.mixed <- any.cont & any.categ

  only.categ <- !any.cont
  only.cont  <- !any.categ


  # Prepare list of variables names for comparison with bic's output names __________________________________________


  factors.preferred.colnames <- gsub(":", "*", vars.in.model)
  
  # Changed in version 1.4.1
  #tmp.expected.varnames  <- gsub(":", ".", vars.in.model)
  tmp.expected.varnames <- vars.in.model


  if (namesx.are.short)
  {
    preferred.varnames <- order1.vars.in.model
    expected.varnames  <- order1.vars.in.model
    preferred.varnames.nolabels <- preferred.varnames
    components <- preferred.varnames


    for (i in which(order>1))
    {
      expected.varname  <- tmp.expected.varnames[i]
      preferred.varname <- factors.preferred.colnames[i]

      if (any.categ[i])
      {
        if (any.cont[i])
        {
          tmp.vars <- names(factors[,i][factors[,i]>0])
          tmp.varnos <- match(tmp.vars, vars.in.model)
          expected.varname  <- character()
          preferred.varname <- character()

          for (j in tmp.varnos)
          {
            if (any.cont[j])
            {
              if (length(expected.varname) == 0)
              {
                expected.varname  <- vars.in.model[j]
                preferred.varname <- vars.in.model[j]
              }
              else
              {
                expected.varname  <- paste(expected.varname,  vars.in.model[j], sep=':')
                preferred.varname <- paste(preferred.varname, vars.in.model[j], sep='*')
              }
            }
            else
            {
              tmp.expected.varname  <- paste(vars.in.model[j], unlist(order1.factors.levels[vars.in.model[j]]), sep='')
              tmp.preferred.varname <- paste(vars.in.model[j], unlist(order1.factors.levels[vars.in.model[j]]), sep='=')

              if (length(expected.varname) == 0)
              {
                expected.varname  <- tmp.expected.varname
                preferred.varname <- tmp.preferred.varname
              }
              else
              {
                expected.varname  <- paste(rep(expected.varname,  rep(length(tmp.expected.varname),  length(expected.varname))),  tmp.expected.varname,  sep=':')
                preferred.varname <- paste(rep(preferred.varname, rep(length(tmp.preferred.varname), length(preferred.varname))), tmp.preferred.varname, sep='*')
              }
            }
          }
        }
        else
        {
          expected.varname <- paste(expected.varname, '..', sep='')
        }
      }

      # add a .1 extension if expected varname already exists in order-1 variables

      tmp <- match(expected.varname, order1.vars.in.model)
      ext <- ifelse(is.na(tmp), "", ".1")
      expected.varname <- paste(expected.varname, ext, sep='')

      expected.varnames  <- c(expected.varnames,  expected.varname)
      preferred.varnames <- c(preferred.varnames, preferred.varname)
      preferred.varnames.nolabels <- c(preferred.varnames.nolabels, rep(factors.preferred.colnames[i], length(preferred.varname)))
    }
  }
  else
  {
    # namesx.are.short = F 
    # (bic.out$factor.type == F || bicreg output)

    my.tmp.l <- tmp.l[order1.vars.in.model]

    expected.varnames  <- c(names(my.tmp.l)[my.tmp.l==0], tmp.longnames)
    preferred.varnames <- c(names(my.tmp.l)[my.tmp.l==0], tmp.preferred.longnames)
    preferred.varnames.nolabels <- c(names(my.tmp.l)[my.tmp.l==0], tmp.names)
    components <- preferred.varnames


    for (i in which(order>1))
    {
      expected.varname  <- tmp.expected.varnames[i]
      preferred.varname <- factors.preferred.colnames[i]

      if (any.categ[i])
      {
        tmp.vars <- names(factors[,i][factors[,i]>0])
        tmp.varnos <- match(tmp.vars, vars.in.model)

        if (any.cont[i])
        {
          expected.varname  <- character()
          preferred.varname <- character()

          for (j in tmp.varnos)
          {
            if (any.cont[j])
            {
              if (length(expected.varname) == 0)
              {
                expected.varname  <- vars.in.model[j]
                preferred.varname <- vars.in.model[j]
              }
              else
              {
                expected.varname  <- paste(expected.varname,  vars.in.model[j], sep='.')
                preferred.varname <- paste(preferred.varname, vars.in.model[j], sep='*')
              }
            }
            else
            {
              tmp.expected.varname  <- paste(vars.in.model[j], unlist(order1.factors.levels[vars.in.model[j]]), sep='')
              tmp.preferred.varname <- paste(vars.in.model[j], unlist(order1.factors.levels[vars.in.model[j]]), sep='=')

              if (length(expected.varname) == 0)
              {
                expected.varname  <- tmp.expected.varname
                preferred.varname <- tmp.preferred.varname
              }
              else
              {
                expected.varname  <- paste(rep(expected.varname,  rep(length(tmp.expected.varname),  length(expected.varname))),  tmp.expected.varname,  sep='.')
                preferred.varname <- paste(rep(preferred.varname, rep(length(tmp.preferred.varname), length(preferred.varname))), tmp.preferred.varname, sep='*')
              }
            }
          }
        }
        else
        {
          Tmp.expected.varname <- paste(expected.varname, '..', sep='')
          expected.varname <- character()

          for (j in tmp.varnos)
          {
            tmp.expected.varname  <- paste(vars.in.model[j], unlist(order1.factors.levels[vars.in.model[j]]), sep='')
            tmp.preferred.varname <- paste(vars.in.model[j], unlist(order1.factors.levels[vars.in.model[j]]), sep='=')

            if (length(expected.varname) == 0)
            {
              expected.varname  <- tmp.expected.varname
              preferred.varname <- tmp.preferred.varname
            }
            else
            {
              expected.varname  <- paste(rep(expected.varname,  rep(length(tmp.expected.varname),  length(expected.varname))),  tmp.expected.varname,  sep=':')
              preferred.varname <- paste(rep(preferred.varname, rep(length(tmp.preferred.varname), length(preferred.varname))), tmp.preferred.varname, sep='*')
            }
          }

          expected.varname <- c(expected.varname, 'ref')
          expected.varname <- paste(Tmp.expected.varname, expected.varname, sep='')

          tmp.ref <- unlist(lapply(order1.factors.levels[tmp.vars],reference.value))
          tmp.ref <- paste(paste(names(tmp.ref), tmp.ref, sep='='), collapse=' or ')
          preferred.varname <- c(preferred.varname, paste(factors.preferred.colnames[i], bic.ref, tmp.ref, ')', sep=''))
        }
      }

      expected.varnames  <- c(expected.varnames,  expected.varname)
      preferred.varnames <- c(preferred.varnames, preferred.varname)
      preferred.varnames.nolabels <- c(preferred.varnames.nolabels, rep(factors.preferred.colnames[i], length(preferred.varname)))
    }
  }


  # complete the list of components (using labels in mixed interaction terms, and in categorical variables if bic.out$factor.type=F)

  for (i in which(order>1))
  {
    tmp.vars <- names(factors[,i][factors[,i]>0])
    tmp.varnos <- match(tmp.vars, vars.in.model)

    components.i <- matrix("", nrow=1, ncol=0)

    tmp.use.labels <- (any(any.categ[tmp.varnos]) & any(any.cont[tmp.varnos])) | (!any(any.cont[tmp.varnos]) & !bic.out.factor.type)

    if (tmp.use.labels)
    {
      for (j in tmp.varnos)
      {
        if (any.categ[j])
        {
          components.j <- paste(vars.in.model[j], unlist(order1.factors.levels[vars.in.model[j]]), sep='=')
        }
        else
        {
          components.j <- vars.in.model[j]
        }

        components.i <- cbind(components.i[rep(seq(nrow(components.i)), rep(length(components.j), nrow(components.i))),], components.j)
      }

      components.i <- lapply(apply(components.i, 1, list), unlist)
    }
    else
    {
      components.i <- list(tmp.vars)
    }
    
    add.ref <- !any(any.cont[tmp.varnos]) & !bic.out.factor.type
    if (add.ref) components.i <- c(components.i, NA)

    components <- c(components, components.i)
  }


  # crash if any of the expected variables names is duplicated AND present in $namesx

  duplicated.varnames <- unique(expected.varnames[duplicated(expected.varnames)])
  duplicated.varnames <- intersect(duplicated.varnames, bic.out.namesx)

  if (length(duplicated.varnames) > 0)
  {
    stop.msg <- "Variables names saved to $namesx were not clear: bic.process() prefers to stop, as it can't know:"

    for (v in duplicated.varnames)
    {
      tmp <- preferred.varnames[which(!is.na(match(expected.varnames, v)))]
      lov <- paste(c('(', paste(tmp, collapse='|'), ')'), collapse='')
      tmp <- paste('which of', lov, 'was assigned the name', v, collapse=' ')
      stop.msg <- c(stop.msg, tmp)
    }

    stop.msg <- paste(stop.msg, collapse='\n')
    stop(stop.msg)
  }

  if (bic.fct == "bicreg")
  {   
    tmp.names  <- names(order1.factors.levels)
    tmp.l      <- pmax(unlist(lapply(order1.factors.levels, length))-1, 1)
    tmp.f      <- order1.factors.levels
    remove.1st <- function(x){x[-1]}
    tmp.f      <- lapply(tmp.f, remove.1st)
    tmp.w      <- which(unlist(lapply(tmp.f,is.null)))
    tmp.f[tmp.w] <- ""
    potential.output.names <- paste(rep(tmp.names, tmp.l), unlist(tmp.f), sep='')
    bicreg.ncateg <- tmp.l
     
    new.output.names <- order1.factors.levels
  }
  else
  {
    new.output.names     <- bic.out$output.names
    names(new.output.names) <- bic.out.namesx
    potential.output.names  <- names(new.output.names)
  }
  
  
  expected.varnames.matchj <- match(potential.output.names, expected.varnames)


  if (any(is.na(expected.varnames.matchj)))
  {
    var <- paste("'", potential.output.names[which(is.na(expected.varnames.matchj))[1]], "'", sep='')
    lov <- paste('(', paste(factors.preferred.colnames, collapse='|'), ')', sep='')

    stop.msg <- paste("Could not match", var, "(from $namesx) with any of the variables in\n  ", lov, sep=' ')
    stop(stop.msg)
  }


  preferred.names <- preferred.varnames[expected.varnames.matchj]
  
  if (bic.fct != "bicreg") names(new.output.names) <- preferred.names

  preferred.names.nolabels <- preferred.varnames.nolabels[expected.varnames.matchj]
  matching.components <- components[expected.varnames.matchj]


  higher.order.only.categ <- order > 1 & !any.cont

  factors.colnos <- match(preferred.names.nolabels, factors.preferred.colnames)
  output.names.only.categ <- higher.order.only.categ[factors.colnos]

  # function will not deal with bic.glm or bic.surv outputs with interactions between categorical variables
  # if factor.type = F if user wants to discard incomplete models

  crash <- discard.incomplete.models & !bic.out.factor.type & any(order[factors.colnos] > 1 & only.categ[factors.colnos]) 

  if (crash) stop("bic.process cannot be run on outputs obtained with discard.incomplete.models=T and factor.type=F when there is one or more interaction term between categorical variables in formula.\nSorry.\n") 


  # change value labels in $output.names (relevant only when bic.out$factor.type==T)

  if (bic.out.factor.type)   
  {
    for (i in which(output.names.only.categ))
    {
      factor.colno <- factors.colnos[i]
      tmp.vars <- names(factors[,factor.colno][factors[,factor.colno]>0])

      expected.labels  <- character()
      preferred.labels <- character()

      for (var in tmp.vars)
      {
        labels <- levels(data[,var])
        tmp.labels <- paste(var, labels, sep='')

        if (length(expected.labels) == 0)
        {
          expected.labels  <- tmp.labels
          preferred.labels <- labels
        }
        else
        {
          expected.labels  <- paste(rep(expected.labels,  rep(length(tmp.labels), length(expected.labels))),  tmp.labels, sep=':')
          preferred.labels <- paste(rep(preferred.labels, rep(length(labels),     length(preferred.labels))), labels,     sep='::')
        }
      }
      
      if (length(tmp.vars) > 1)
      {
        # change 'ref' into the preferred long label
        
        expected.labels <- c(expected.labels, 'ref')
        
        tmp.ref <- unlist(lapply(order1.factors.levels[tmp.vars],reference.value))
        tmp.ref <- paste(paste(names(tmp.ref), tmp.ref, sep='='), collapse=' or ')
        tmp.ref <- paste(substr(bic.ref,2,nchar(bic.ref)), tmp.ref, ')', sep='')
        preferred.labels <- c(preferred.labels, tmp.ref)
      }

      tmp <- match(new.output.names[[i]], expected.labels)
      w <- which(!is.na(tmp))
      new.output.names[[i]][w] <- preferred.labels[tmp[w]]
    }
  }

  # find which vars are categorical

  categ.vars.order1 <- vars.in.model[any.categ & order==1]
  if (length(categ.vars.order1)) extra$reference.values <- unlist(lapply(order1.factors.levels[categ.vars.order1], reference.value), use.names=T)


  # Define mle colnames ___________________________________________________________________________


  if (namesx.are.short)
  {
    extract.output.names <- function(x) ifelse(is.na(x[2]),x[1],paste(x,collapse='='))

    labels <- unlist(new.output.names, use.names=F)
    l <- unlist(lapply(new.output.names, length), use.names=F)
    names <- rep(preferred.names, l)
    tmp.names <- rep(preferred.names.nolabels, l)
    keep <- duplicated(names) | is.na(labels)

    mle.names0 <- names[keep]
    mle.names0.nolabels <- tmp.names[keep]
    tmp <- matrix(c(mle.names0,labels[keep]), ncol=2)
    mle.names <- apply(tmp, 1, extract.output.names)
  }
  else
  {
    mle.names <- preferred.names
    mle.names0 <- preferred.names
    mle.names0.nolabels <- preferred.names.nolabels
  }


  if (bic.fct == "surv") intercept.in.model <- F


  if (intercept.in.model) 
  {
    mle.names  <- c("Intercept", mle.names)
    mle.names0 <- c(NA, mle.names0)
    mle.names0.nolabels <- c(NA, mle.names0.nolabels)
  }

  # ______________________________________________________________________________________________________________

  which <- bic.out$which 
  
  # ______________________________________________________________________________________________________________

  # bicreg doesn't seem to eliminate incomplete models with regards to factor variables!
  # We will do that now.
  
  if (bic.fct == "bicreg" && any(only.categ) && bic.out.factor.type)
  {
    tmp <- matrix(categ.vars.order1, byrow=T, ncol=length(categ.vars.order1), nrow=ncol(which))
    tmp <- preferred.names.nolabels == tmp
    tmp <- t(which %*% tmp)
    tmp <- tmp==0 | tmp==bicreg.ncateg[categ.vars.order1]
    order1.ok <- apply(tmp,2,all)
    
    which <- which[,!duplicated(preferred.names.nolabels), drop=F]
    colnames(which) <- unique(preferred.names.nolabels)
  }
  else
  {
    order1.ok <- rep(T, nrow(which))
    colnames(which) <- preferred.names
  }

  # ______________________________________________________________________________________________________________
    
  mixed.interactionterms.factor.type <- mixed.interactionterms.factor.type & bic.out.factor.type
    # (forced to be false when bic.out$factor.type was F)


  if (mixed.interactionterms.factor.type && any(any.mixed))
  {
    # discard models that are incomplete due to absence of categories in mixed interaction terms 
    # (regardless of absence/presence of lower-order terms at this point)

    tmp.names <- factors.preferred.colnames[any.mixed]
    tmp.col <- match(preferred.names.nolabels, tmp.names)
    tmp.row <- which(!is.na(tmp.col))
    tmp.col <- tmp.col[tmp.row]

    tmp.which <- matrix(F, nrow=ncol(which), ncol=length(tmp.names))
    tmp.which[(tmp.col-1)*nrow(tmp.which) + tmp.row] <- T
    tmp.nterms <- as.vector(table(tmp.col))

    mixedterms.ok <- t(which %*% tmp.which)
    mixedterms.ok <- mixedterms.ok == 0 | mixedterms.ok == tmp.nterms
    mixedterms.ok <- apply(mixedterms.ok, 2, all)
  }
  else
  {
    mixedterms.ok <- rep(T, nrow(which))
  }

  # ______________________________________________________________________________________________________________


  if (discard.incomplete.models && any(order>1))
  {
    # if any higher-order term (that is, any interaction term) is present in formula, 
    # we may need to clean/reduce the list of models

    tmp.m <- ncol(which)

    child.mixed  <- matrix(any.mixed[factors.colnos], nrow=tmp.m, ncol=tmp.m)
    parent.mixed <- t(child.mixed)

    child.only.categ  <- matrix(only.categ[factors.colnos], nrow=tmp.m, ncol=tmp.m)
    parent.only.categ <- t(child.only.categ)

    tmp.parent.cols <- which(order[factors.colnos] > 1)

    parent.mixed      <- parent.mixed[     , tmp.parent.cols, drop=F]
    parent.only.categ <- parent.only.categ[, tmp.parent.cols, drop=F]
    child.mixed       <- child.mixed[      , tmp.parent.cols, drop=F]
    child.only.categ  <- child.only.categ[ , tmp.parent.cols, drop=F]

    if (bic.out.factor.type)
    {
      look.at.components.match <- parent.mixed & child.mixed
    }
    else
    {
      look.at.components.match <- parent.mixed | child.only.categ
    }

    # components.match
   
    components.match <- matrix(T, nrow=tmp.m, ncol=0)

    for (k in tmp.parent.cols)
    {
      tmp.components.match <- unlist(lapply(lapply(matching.components, match,table=matching.components[[k]],nomatch=0),min))>0
      components.match <- cbind(components.match, tmp.components.match)
    }

    # names.match

    tmp.nvars <- apply(factors, 2, sum)
    tmp <- t(factors) %*% factors
    tmp <- tmp == tmp.nvars
    names.match <- tmp[, factors.colnos, drop=F]
    names.match <- names.match[factors.colnos,, drop=F]
    names.match <- names.match[, order[factors.colnos] > 1, drop=F]

    # final match

    lower.order.terms.included <- (look.at.components.match & components.match) | (!look.at.components.match & names.match)
    # indicates whether ALL lower order terms AND interaction terms are present, model by model    

    lower.order.terms.m <- apply(lower.order.terms.included, 2, sum)
   
    tmp <- t(which %*% lower.order.terms.included)

    tmp.absent <- diag(tmp.m)[,tmp.parent.cols,drop=F]
    tmp.absent <- t(which %*% tmp.absent)    

    ok.interns <- (tmp==lower.order.terms.m) | (tmp.absent == 0)
    ok.interns <- apply(ok.interns, 2, all)
  }
  else
  {
    ok.interns <- rep(T, nrow(which))
  }

  # _____________________________________________________________________


  which.ok <- mixedterms.ok & ok.interns & order1.ok

  exclusion.criterion <- pmax(!mixedterms.ok, 2*!ok.interns, 3*!order1.ok)
  exclusion.criterion <- as.character(factor(exclusion.criterion,levels=seq(3),labels=c('dropped level in mixed interaction term','interaction term without lower order terms being in model','dropped factor level')))

  models2keep <- which(which.ok)

  if (length(models2keep)==0) 
    stop("No model left after discarding incomplete/incorrect models.")
    
  exclusion.criterion <- exclusion.criterion[!is.na(exclusion.criterion)]
  
  if (length(exclusion.criterion) > 0)
  {
    excluded <- which(!which.ok)
    which.excluded <- which[!which.ok,,drop=F]
  }
  else
  {
    excluded <- NULL
    which.excluded <- NULL
    exclusion.criterion <- NULL
  }

  OR.corrected <- NULL
  
  if (!which.ok[1])
  {
    tmp <- bic.out$postprob[c(1, models2keep[1])]
    tmp <- tmp/(1-tmp)
    OR.corrected <- tmp[1]/tmp[2] * OR
  }

  changes <- any(!which.ok)

  if (!changes)
  {
    extra$nests <- bic.out$nests
    extra$disp  <- bic.out$disp
  }
  # else we ignore these objects, which would (probably) need recalculation but we find not interesting...


  which          <- which[models2keep,,drop=F]
  bic            <- bic.out$bic[models2keep]
  mle            <- bic.out$mle[models2keep,,drop=F]
  se             <- bic.out$se[models2keep,,drop=F]
  postprob       <- bic.out$postprob[models2keep]
  n.models       <- length(models2keep)
  
  if (bic.fct == "bicreg") r2 <- bic.out$r2[models2keep]

  size                <- bic.out$size[models2keep]
  prior.model.weights <- bic.out$prior.model.weights[models2keep]
  deviance            <- bic.out$deviance[models2keep]


  if (bic.out.factor.type)
  {
    if (bic.fct == "bicreg")
    {
      mle.whichCol <- match(mle.names0.nolabels, colnames(which))
    }
    else
    {
      mle.whichCol <- match(mle.names0, colnames(which))
    }
  }
  else
  {
    mle.whichCol <- match(mle.names0, preferred.names)
  }
  
  mle.which <- which[,mle.whichCol,drop=F]
  if (intercept.in.model) mle.which[,1] <- T


  if (changes)
  {
    postprob <- postprob/sum(postprob)
    postprob.as.row <- matrix(postprob, nrow=1)

    probne0  <- as.vector(postprob.as.row %*% which)

    postmean   <- as.vector(postprob.as.row %*% mle)
    postsqmean <- as.vector(postprob.as.row %*% (mle^2 + se^2))
    postsd     <- sqrt(postsqmean-postmean^2)

    # cond means and se's

    cond.probne0 <- as.vector(postprob.as.row %*% mle.which)
    if (intercept.in.model) cond.probne0[1] <- 1

    condpostmean   <- as.vector(postprob.as.row %*% (mle.which * mle))           / cond.probne0
    condpostsqmean <- as.vector(postprob.as.row %*% (mle.which *(mle^2 + se^2))) / cond.probne0

    condpostsd     <- sqrt(condpostsqmean-condpostmean^2)

    probne0  <- round(100*probne0, digits=1)
  }
  else
  {
    probne0 <- bic.out$probne0

    postmean       <- bic.out$postmean
    postsd         <- bic.out$postsd
    condpostmean   <- bic.out$condpostmean
    condpostsd     <- bic.out$condpostsd
  }


  # --------------------------------------------------------------------------------------
  # See for which terms we do not compute exp(mle) 
  # [either for odds ratios, rate ratios or hazard ratios]
  # we do not compute exp(mle|mle.lower|mle.upper) for variable x when it is
  # involved in a higher-order interaction term


  compute.exp <- matrix(T, nrow=nrow(mle), ncol=ncol(mle)-intercept.in.model)

  if (any(order[factors.colnos] > 1))
  {
    tmp.w  <- which(order[factors.colnos]>1)
    compute.exp[,tmp.w] <- F
    
    donot.compute.exp <- (which[,tmp.w,drop=F] %*% t(factors[,factors.colnos[tmp.w],drop=F])) > 0
    tmp.w <- which(donot.compute.exp)
    tmp.r <- row(donot.compute.exp)[tmp.w]
    tmp.c <- col(donot.compute.exp)[tmp.w]
    tmp.m <- match(colnames(donot.compute.exp),preferred.names.nolabels)
    tmp.c <- tmp.m[tmp.c]
    tmp.index <- (tmp.c-1)*nrow(which) + tmp.r
    compute.exp[tmp.index] <- F
  }

  # --------------------------------------------------------------------------------------

  colnames(mle) <- mle.names
  colnames(se)  <- mle.names

  if (recode.absent.var.value.as.na)
  {
    # put NA's in mle's and se's where variables are not in model (if indicated)

    parm.not.in.model <- which(!mle.which)

    mle[parm.not.in.model] <- NA
    se[parm.not.in.model]  <- NA  
  }

  if (report.ci)
  {  
    z <- qnorm((1+ci.level)/2)
    mle.lower <- mle - z*se
    mle.upper <- mle + z*se
  }
  else
  {
    mle.lower <- NULL
    mle.upper <- NULL
  }


  odds.ratio       <- NULL
  odds.ratio.lower <- NULL
  odds.ratio.upper <- NULL

  if (length(report.or) && report.or)
  {
    odds.ratio <- exp(mle)
    if (intercept.in.model) odds.ratio <- odds.ratio[,-1,drop=F]
    odds.ratio[!compute.exp] <- NA

    if (report.ci)
    {
      odds.ratio.lower <- exp(mle.lower)
      if (intercept.in.model) odds.ratio.lower <- odds.ratio.lower[,-1,drop=F]
      odds.ratio.lower[!compute.exp] <- NA

      odds.ratio.upper <- exp(mle.upper)
      if (intercept.in.model) odds.ratio.upper <- odds.ratio.upper[,-1,drop=F]
      odds.ratio.upper[!compute.exp] <- NA
    }
  }


  hazard.ratio <- NULL
  hazard.ratio.lower <- NULL
  hazard.ratio.upper <- NULL

  if (length(report.hazard.ratio) && report.hazard.ratio)
  {
    hazard.ratio <- exp(mle)
    hazard.ratio[!compute.exp] <- NA

    if (report.ci)
    {
      hazard.ratio.lower <- exp(mle.lower)
      hazard.ratio.lower[!compute.exp] <- NA
      
      hazard.ratio.upper <- exp(mle.upper)
      hazard.ratio.upper[!compute.exp] <- NA
    }
  }

  rate.ratio <- NULL
  rate.ratio.lower <- NULL
  rate.ratio.upper <- NULL

  if (length(report.rate.ratio) && report.rate.ratio)
  {
    rate.ratio <- exp(mle)
    rate.ratio[!compute.exp] <- NA

    if (report.ci)
    {
      rate.ratio.lower <- exp(mle.lower)
      rate.ratio.lower[!compute.exp] <- NA
      
      rate.ratio.upper <- exp(mle.upper)
      rate.ratio.upper[!compute.exp] <- NA
    }
  }

  # Attach names to some objects

  if (bic.fct != "bicreg")
  {
    prior.param <- bic.out$prior.param
    names(prior.param) <- preferred.names
  }
  else
  {
    prior.param <- NULL
  }

  names(probne0)      <- colnames(which)

  names(postmean)     <- mle.names
  names(postsd)       <- mle.names
  names(condpostmean) <- mle.names
  names(condpostsd)   <- mle.names


  # --------------------------------------------------------------------------------------------------------------
  # If formula included interaction terms involving continuous variables, factor.type was T and
  # incomplete models were discarded, then a few objects need to be changed.


  simplify.internterms.with.contvars <- discard.incomplete.models & any(order>1 & any.cont) & bic.out.factor.type

  if (simplify.internterms.with.contvars)
  {
    tmp.o <- cumsum(!duplicated(preferred.names.nolabels))
    tmp.drop <- duplicated(preferred.names.nolabels) & order[tmp.o] > 1 & any.cont[tmp.o]

    tmp.colnames <- preferred.names
    tmp.w <- which(tmp.drop)
    tmp.w <- which(!is.na(match(preferred.names.nolabels, preferred.names.nolabels[tmp.w])))
    tmp.colnames[tmp.w] <- preferred.names.nolabels[tmp.w]
    
    colnames(which) <- tmp.colnames
    which <- which[, !tmp.drop, drop=F]

    names(prior.param) <- tmp.colnames
    prior.param <- prior.param[!tmp.drop]

    names(probne0) <- tmp.colnames
    probne0 <- probne0[!tmp.drop]
  }


  first.model.in  <- first.row(which)
  first.model.out <- first.row(!which)
  
  discarded.vars <- which(is.na(first.model.in))
  preferred.vars.in.model <- gsub(":", "*", vars.in.model)

  # Redefine models labels

  model.label <- function(w){paste(names(which(w)),collapse=",")}
  label <- apply(which, 1, model.label)
  label[nchar(label)==0] <- "NULL"


  # --- Join results in a unique list, where undefined objects are not reported -----------------------------


  out <- list(label = label, bic = bic, n.models = n.models, probne0 = probne0,
    which = which, size = size, mle = mle, se = se,
    postmean = postmean, postsd = postsd,
    postprob = postprob,
    condpostmean = condpostmean, condpostsd = condpostsd,
    n = n, vars.in.model = preferred.vars.in.model, first.model.in = first.model.in, first.model.out = first.model.out,
    input.nrows = nrow(data), 
    reduced = bic.out$reduced,
    call = bic.out$call, factor.type = bic.out.factor.type,
    original.namesx = bic.out$namesx,
    discard.incomplete.models = discard.incomplete.models, mixed.interactionterms.factor.type = mixed.interactionterms.factor.type
    )


  tmp <- preferred.vars.in.model[prior.param==1]
  if (length(tmp)) extra$forced.in <- tmp
  
  extra$family              <- bic.out$family
  extra$dropped             <- bic.out$dropped
  extra$OR.corrected        <- OR.corrected
  extra$prior.model.weights <- prior.model.weights
  extra$prior.param         <- prior.param

  extra$mle.lower <- mle.lower
  extra$mle.upper <- mle.upper

  extra$odds.ratio       <- odds.ratio
  extra$odds.ratio.lower <- odds.ratio.lower
  extra$odds.ratio.upper <- odds.ratio.upper

  extra$hazard.ratio       <- hazard.ratio
  extra$hazard.ratio.lower <- hazard.ratio.lower
  extra$hazard.ratio.upper <- hazard.ratio.upper

  extra$rate.ratio       <- rate.ratio
  extra$rate.ratio.lower <- rate.ratio.lower
  extra$rate.ratio.upper <- rate.ratio.upper

  extra$censoring.var    <- censoring.var
  extra$rejected.vars    <- bic.out$rejected.vars
  extra$notes            <- notes

  if (length(n.events) > 0) extra$n.events <- n.events
  
  if (report.r2)
  {
    if (bic.fct == "bicreg")
    {
      extra$r2 <- r2/100
    }
    else
    {
      SST <- sum((bic.out$y-mean(bic.out$y))^2)
      extra$r2 <- 1 - deviance/SST
    }
  }
  
  extra$deviance <- deviance

  extra$excluded            <- excluded
  extra$which.excluded      <- which.excluded
  extra$exclusion.criterion <- exclusion.criterion

  append(out, extra)
} # end of bic.process
