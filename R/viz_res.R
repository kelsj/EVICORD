#' function: viz_res.R
#' takes list of output from gibbs sampler, visualize posterior chains
#' input: filename e.g. gibbs_sampler_out.rds, ac (allele count), tf (integer factor to thin chains by), plot width (optional, default 4in), plot height (optional, default 8in)
#' requires ggplot, cowplot, dplyr

#' output plots:
#' 1. chain of beta samples by iteration
#' 2. chain of pi samples by iteration
#' 3. density of t estimates for each allele pair category (t1,t2,t3)
#' 4. empirical cumulative density of the fraction of z samples for each variant greater than x (with x=1,2,3,...); i.e., zGT1 = the fraction of posterior samples that the variant is assigned to be non-IBD

#function to get summary of sampler chains from rds file
viz_res = function(filename, ac, tf, width=4, height=8, outfile="gibbs_sampler_out_viz.pdf"){
	npairs = choose(ac,2)
	nk = as.integer(ac/2)+1

	#data = list("nVars"=n,"beta"=beta_i,"pi"=pi_chain,"t1"=t1_i,"t2"=t2_i,"t3"=t3_i,"z"=z_chain,"v"=v_chain)
	output = readRDS(filename)
	
	nvars = output$nVars
	
	beta = output$beta %>% data.frame() %>% tbl_df()
	niter=nrow(beta)-1
	beta$iter = 1:(niter+1)
	
	thin = seq(2,(niter+1)/tf,1)*tf
	medBeta = median(beta[thin,]$X1)
	
	p1 = ggplot(beta %>% filter(iter %in% thin)) + geom_line(aes(x=iter,y=X1)) + ylab("IBD beta")
 	
	pi=output$pi%>% data.frame() %>% tbl_df()
	colnames(pi) = paste("k",1:nk,sep="")
	pi$iter = 1:(niter+1)
	pi = pi %>% gather(key="k",value="pi",k1:paste("k",nk,sep=""))
	p2 = ggplot(pi %>% filter(iter %in% thin)) + geom_line(aes(x=iter,y=pi,color=k))
	
	t1=output$t1%>% data.frame() %>% tbl_df()
	t1$iter = 1:(niter+1)
	lastcol = paste("X",nvars,sep="")
	t1 = t1 %>% gather(key="var",value="tEst",X1:!!lastcol) %>% mutate(t="t1")
	t2=output$t2%>% data.frame() %>% tbl_df()
	t2$iter = 1:(niter+1)
	t2 = t2 %>% gather(key="var",value="tEst",X1:!!lastcol) %>% mutate(t="t2")
	t3=output$t3%>% data.frame() %>% tbl_df()
	t3$iter = 1:(niter+1)
	t3 = t3 %>% gather(key="var",value="tEst",X1:!!lastcol) %>% mutate(t="t3")
	t = bind_rows(t1,t2,t3)
	tsummary = t %>% filter(iter %in% thin) %>% group_by(var,t) %>% summarise(tMed = median(tEst,na.rm=T),mean=mean(tEst,na.rm=T))
	p3 = ggplot(tsummary) + geom_density(aes(x=tMed,color=t)) + scale_x_continuous(limits=c(0,1000))
	
	z = output$z %>% data.frame() %>% tbl_df()

	fgt = c(1:nvars)
	for (x in 1:(nk-1)){
		zgtx = apply(z[thin,],2,fracKgtX,x=x)
		fgt = cbind(fgt,zgtx)
	}	
	colnames(fgt) = c("varNum",paste("zGT",1:(nk-1),sep=""))
	fgt = tbl_df(fgt)
	fgt = fgt %>% gather(key="k",value="zGT",zGT1:!!paste("zGT",(nk-1),sep=""))
	p4 = ggplot(fgt) + stat_ecdf(aes(x=zGT,color=k))
	plotsummary = plot_grid(p1,p2,p3,p4,nrow=4)	
	ggsave(outfile,plotsummary,width=width,height=height)
}

