#! /usr/bin/env Rscript

Nonpareil.f <- function(
    ### Function of the projected model.
    x,
    ### Values of sequencing effort (typically in bp).
    a,
    ### Parameter alpha of the gamma CDF.
    b
    ### Parameter beta of the Gamma CDF.
    ){ return(pgamma(log(x+1),a, b))
    ### Predicted values of abundance-weighted average coverage.
    }
Nonpareil.antif <- function(
    ### Complement function of Nonpareil.f.
    y,
    ### Values of abundance-weighted average coverage
    a,
    ### Parameter alpha of the gamma CDF.
    b
    ### Parameter beta of the gamma CDF.
    ){ return(exp(qgamma(y,a,b))-1)
    ### Estimated sequencing effort.
    }

simple_curve <- function(
    ### Generates a Nonpareil curve from a .npo file.
    file,
    factor=1,
    star=95,
    plot=TRUE,
    ### 
    r=NA,
    ### Red component of the curve's color. If NA, it's randomly set. If <=1, it's assumed to be in the range [0,1]; if
    ### >1, it's assumed to be in the range [0,256].
    g=NA,
    ### Green component of the curve's color. Same as `r`.
    b=NA,
    ### Blue component of the curve's color. Same as `r`.
    curve.lwd=2,
    ### Line width of the Nonpareil curve.
    curve.alpha=0.4,
    ### Alpha value (from 0 to 1) of the Nonpareil curve.
    model.lwd=1,
    ### Line width of the model.
    model.alpha=1,
    ### Alpha value (from 0 to 1) of the model.
    log='x',
    ### Axis to plot in logarithmic scale. It can be 'x' (sequencing effort, default), 'y' (coverage), 'xy' (both logarithmic),
    xmax=10e12,
    ### Maximum sequencing effort to plot.
    xmin=1e3,
    ### Minimum sequencing effort to plot
    ymax=1,
    ### Maximum coverage to plot.
    ymin=1e-6,
    logsampling=0,
    weights.exp=NA
    ### Vector of values to be tested (in order) as exponent of the weights distribution. If the model fails to converge, sometimes
    ### manual modifications in this parameter may help. By default (NA), five different values are tested in the following order:
    ### For linear sampling, -1.1, -1.2, -0.9, -1.3, -1. For logarithmic sampling (-d option in Nonpareil), 0, 1, -1, 1.3, -1.1.
    ){
    # Examine consistency
    if(is.na(file)) stop('The .npo file argument is mandatory');

    # Read input
    out <- list(kappa=0, C=0, LRstar=0, LR=0, modelR=0, diversity=0);
    a <- read.table(file, sep="\t", h=F);
    out$kappa <- a$V2[nrow(a)];
    a$V1 = a$V1*factor;
    num_reads = max(a$V1);
    LR <- num_reads;   # could also number total bp
    out$LR <- LR;
    values = a$V2; # take mean; options are mean,q1,q2,q3
    out$C <- max(values);
    twenty.pc = which.max(a$V1[a$V1<=0.2*max(a$V1)]);
    if(a[twenty.pc, 5]==0){
        warning('**Median of the curve is zero at 20% of the reads,\
                sequencing depth is too low.');
        }
    if(a[nrow(a), 2]<=1e-5){
        warning('**Curve too low. Sequencing depth is too low.');
        }
    if(a[2, 2]>=1-1e-5){
        warning('Curve too steep, check parameters and re-run \
                (e.g., decrease sumsampling interval).');
        }

    # Draw it
    if(plot){
        figfile<-paste(file, '.R.pdf', sep='');
        pdf(figfile);
        ylab <- 'Estimated average coverage';
        xlab <- 'Number of reads collected';
        plot(1, t='n', xlab=xlab, ylab=ylab, log=log, xlim=c(xmin, xmax), ylim=c(ymin, ymax));
        abline(h=1, lty=2, col='gray80');
        abline(h=star/100, lty=2, col='green');
        abline(v=10^seq(0,15,by=3), lty=2, col='gray80')

        if(is.na(r)) r <- sample(200,1);
        if(is.na(g)) g <- sample(200,1);
        if(is.na(b)) b <- sample(200,1);
        if(r>1) r <- r/256;
        if(g>1) g <- g/256;
        if(b>1) b <- b/256;
        
        err.y <- c(values+a$V3, rev(values-a$V3));
        polygon(c(a$V1, rev(a$V1)), ifelse(err.y<=ymin*0.1, ymin*0.1, err.y), col=rgb(r,g,b,.2), border=NA);
        #polygon(c(a$V1, rev(a$V1)), err.y, col=rgb(r,g,b,.2), border=NA);
        lines(a$V1, values, col=rgb(r,g,b,curve.alpha), lwd=curve.lwd);
        }
    
    # Model it
    sel  <- values>0 & values<0.9;
    x <- a$V1[sel];
    if(length(x)>10){
        y <- values[sel];
        data <- list(x=x, y=y);
        if(is.na(weights.exp[1])){
            # linear sampling
            if(logsampling==0){ weights.exp <- c(-1.1,-1.2,-0.9,-1.3,-1) }
            # log sampling
            else{ weights.exp <- c(0, 1, -1, 1.3, -1.1) }
            }
        weights.i <- 0;
        estModel <- TRUE;
        while(estModel){
            weights.i <- weights.i+1;
            model <- nls(y ~ Nonpareil.f(x, a, b), data=data, weights=(a$V3[sel]^weights.exp[weights.i]), start=list(a=1, b=0.1), lower=c(a=0, b=0), algorithm='port', control=nls.control(minFactor=1e-25000, tol=1e-15, maxiter=1024, warnOnly=T));
            if(is.na(weights.exp[weights.i+1]) | summary(model)$convInfo$isConv) estModel <- FALSE;
            }
        if(summary(model)$convInfo$isConv){
            model.x <- exp(seq(log(xmin), log(xmax), length.out=1e3));
            model.y <- predict(model, list(x=model.x));
            model.x <- model.x;
            if (plot){
                lines(model.x, model.y, col=rgb(r,g,b,model.alpha));
                current_cov = predict(model, list(x=LR))
                points(LR, current_cov, col=rgb(r,g,b), pch=21, bg='white');
                }

            model.x = append(model.x, LR);
            model.y = append(model.y, current_cov);
            df_out = data.frame(model.x, model.y);

            data_outfile<-paste(file, '.model', sep='');
            write.table(df_out, file=data_outfile, sep='\t',
                        row.names=F, col.names=F, quote=F);
            mes <- paste('model data saved in', data_outfile, sep=' ');
            write(mes, stderr())

            out$model.x <- model.x;
            out$model.y <- model.y;
            out$model=model;
            pa <- coef(model)[1];
            pb <- coef(model)[2];
            if(pa > 1) out$diversity <- (pa-1)/pb;
            out$LRstar <- Nonpareil.antif(star/100, pa, pb);
            out$modelR <- cor(y, predict(model, list(x=x)));
            }
        else{
            warning('Model didn\'t converge.');
            }
        }
    else{
        warning('Insufficient resolution below 90% coverage, skipping model');
        }
        mes <- paste('figure saved in', figfile, sep=' ');
        write(mes, stderr())
        dev.off();
    return(out);

    ### A list with the following entries:
    ### 
    ### kappa: "Redundancy" value of the entire dataset.
    ### 
    ### C: Average coverage of the entire dataset.
    ### 
    ### LRstar: Estimated sequencing effort required to reach the objective average coverage (star, 95% by default).
    ### 
    ### LR: Actual sequencing effort of the dataset.
    ### 
    ### modelR: Pearson's R coefficient betweeen the rarefied data and the projected model.
    ###
    ### diversity: Nonpareil-based diversity index. This value's units are the natural logarithm of the units of
    ###    sequencing effort (by default, log-bp), and indicates the inflection point of the fitted model for the
    ###    Nonpareil curve. If the fit doesn't converge, or the model is not estimated, the value is zero (0).
    ### 
    ### model.x (if retturnModelValues=TRUE): Values of sequencing effort in the projected model.
    ### 
    ### model.y (if retturnModelValues=TRUE): Values of average coverage in the projected model.
    ### 
    ### model (if returnModelParameters=TRUE): Fitted model.
    ### 
    }

#debug(simple_curve)
args <- commandArgs(TRUE)
if (length(args) != 2){
    stop('Usage: Rscript <thisFile> <file.npo> factor')
    }

npoFile <- args[1]
factor <- as.numeric(args[2])

out<- simple_curve(npoFile, factor=factor, xmin=1e3, xmax=1e10, logsampling=1)
out_keys<-c('kappa', 'C', 'LRstar', 'LR', 'modelR', 'diversity')
mes = paste(out_keys, out[out_keys], sep='\t', collapse='\n')
#out <- list(kappa=0, C=0, LRstar=0, LR=0, modelR=0, diversity=0);
print(out$model)
write(mes, stdout())
