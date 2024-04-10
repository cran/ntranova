#' Neutrosophic Two Way ANOVA
#'
#' @param dt is a data frame
#'
#' @return Neutrosophic ANOVA Table
#' @export
#' @importFrom stats qf
#' @importFrom utils capture.output
#' @examples
#' y=c(4,5,3,9,11,8,15,12,14)
#' y1=c(6,7,5,11,14,10,17,13,16)
#' tr=c(1,1,1,2,2,2,3,3,3)
#' cek=c(1,2,3,1,2,3,1,2,3)
#' dt=data.frame(y,y1,tr,cek)
#' ntaov(dt)
ntaov=function(dt)
{
  rw=nrow(dt)
  if (is.null(rw) || is.na(rw) || is.infinite(rw)) {
    stop(terminatepg)
  }
  terminatepg = "Terminating program ... \n"
  origmsg = "Here is the original message: "
  rw = tryCatch(rw)
  {
    as.integer(rw)
  }
  warning = function(warning_msg) {
    message("rw The  value is non-numeric, please try again")
    message(origmsg)
    message(paste(warning_msg, "\n"))
    return(NULL)
  }
  dt = tryCatch(dt)
  {
    as.data.frame(dt)
  }
  warning = function(warning_msg) {
    message(" The data frame is empty, please try again")
    message(origmsg)
    message(paste(warning_msg, "\n"))
    return(NULL)
  }

  g=matrix(nrow=2)
  sso=matrix(nrow=2)
  sst=matrix(nrow=2)
  ssr=matrix(nrow=2)
  sse=matrix(nrow=2)
  mst=matrix(nrow=2)
  msr=matrix(nrow=2)
  mse=matrix(nrow=2)
  cf=matrix(nrow=2)
  f1=matrix(nrow=2)
  f2=matrix(nrow=2)
  g[1]=sum(dt[,1])
  n=max(dt[,3])
  m=max(dt[,4])
  nm=nrow(dt)
  str=matrix(nrow=n)
  scek=matrix(nrow=m)
  message("The Basic data","\n")
  message("---------------------------------------------","\n")
  data=matrix(c(dt[,1],dt[,2],dt[,3],dt[,4]),ncol=4)
  colnames(data)=c("Lower bounds","Upper bounds","Treatments","Sectors")
  rownames(data)=c(1:nm)
  message(paste0(capture.output(data),collapse="\n"))
  message("---------------------------------------------","\n")
  for(i in 1:n)
  {
    s=0
    for(j in 1:nm)
    {
      if(dt[j,3]==i)
      {s=s+dt[j,1]}
      str[i]=s
    }
  }
  for(i in 1:m)
  {
    s=0
    for(j in 1:nm)
    {
      if(dt[j,4]==i)
      {
        s=s+dt[j,1]
      }
      scek[i]=s
    }
  }
  cf[1]=g[1]^2/nm
  sso[1]=sum(dt[,1]^2)-cf[1]
  sst[1]=sum(str^2)/(nm/n)-cf[1]
  ssr[1]=sum(scek^2)/(nm/m)-cf[1]
  sse[1]=sso[1]-sst[1]-ssr[1]
  mst[1]=sst[1]/(n-1)
  msr[1]=ssr[1]/(m-1)
  mse[1]=sse[1]/(nm-n-m+1)
  f1[1]=mst[1]/mse[1]
  f1[2]=msr[1]/mse[1]
  g[2]=sum(dt[,2])
  cf[2]=g[2]^2/nm
  for(i in 1:n)
  {
    s=0
    for(j in 1:nm)
    {
      if(dt[j,3]==i)
      {s=s+dt[j,2]}
      str[i]=s
    }
  }
  for(i in 1:m)
  {
    s=0
    for(j in 1:nm)
    {
      if(dt[j,4]==i)
      {
        s=s+dt[j,2]
      }
      scek[i]=s
    }
  }
  sso[2]=sum(dt[,2]^2)-cf[2]
  sst[2]=sum(str^2)/(nm/n)-cf[2]
  ssr[2]=sum(scek^2)/(nm/m)-cf[2]
  sse[2]=sso[2]-sst[2]-ssr[2]
  mst[2]=sst[2]/(n-1)
  msr[2]=ssr[2]/(m-1)
  mse[2]=sse[2]/(nm-n-m+1)
  f2[1]=mst[2]/mse[2]
  f2[2]=msr[2]/mse[2]
  f1=round(f1,2)
  f2=round(f2,2)
  sso=round(sso,2)
  sst=round(sst,2)
  ssr=round(ssr,2)
  sse=round(sse,2)
  mst=round(mst,2)
  msr=round(msr,2)
  mse=round(mse,2)
  if(sst[1]>sst[2])
  {k=sst[1]
  sst[1]=sst[2]
  sst[2]=k}
  if(ssr[1]>ssr[2])
  {k=ssr[1]
  ssr[1]=ssr[2]
  ssr[2]=k}
  if(sso[1]>sso[2])
  {k=sso[1]
  sso[1]=sso[2]
  sso[2]=k}
  if(sse[1]>sse[2])
  {k=sse[1]
  sse[1]=sse[2]
  sse[2]=k}

  if(mst[1]>mst[2])
  {k=mst[1]
  mst[1]=mst[2]
  mst[2]=k}
  if(msr[1]>msr[2])
  {k=msr[1]
  msr[1]=msr[2]
  msr[2]=k}

  if(mse[1]>mse[2])
  {k=mse[1]
  mse[1]=mse[2]
  mse[2]=k}
  if(f1[1]>f2[1])
  {k=f1[1]
  f1[1]=f2[1]
  f2[1]=k}
  if(f1[2]>f2[2])
  {k=f1[2]
  f1[2]=f2[2]
  f2[2]=k}
  #cat(" S.O.V     ","  DF ","           Sum Sq     ","      Men Sq    ","      F    ","\n")
  #cat("-------------------------------------------------------------","\n")
  #cat(" between G.   ",n-1," ","[",sst[1],sst[2],"]" ,"[",mst[1],mst[2],"]","[",f1[1],f2[1],"]","   ","\n")
  #cat(" between cek. ",m-1," ","[",ssr[1],ssr[2],"]" ,"[",msr[1],msr[2],"]","[",f1[2],f2[2],"]","   ","\n")
  #cat(" Residuals    ",nm-n-m+1,"","[",sse[1],sse[2],"]" ,"[",mse[1],mse[2],"]","\n")
  #cat(" Total        ",nm-1,"","[",sso[1],sso[2],"]",  "\n")

  res=matrix(nrow=4,ncol=8)
  rownames(res)=c("Betweem groups","between sectors","Residuals","Total")
  colnames(res)=c("Df","Sum Sq-L","Sum Sq-H","Meam Sq-L","Meam Sq-H","F-L","F-H","Critical F-Value")
  res[1,1]=n-1;res[2,1]=m-1;res[3,1]=nm-n-m+1;res[4,1]=nm-1
  res[1,2]=sst[1];res[2,2]=ssr[1];res[3,2]=sse[1];res[4,2]=sso[1]
  res[1,3]=sst[2];res[2,3]=ssr[2];res[3,3]=sse[2];res[4,3]=sso[2]
  res[1,4]=mst[1];res[2,4]=msr[1];res[3,4]=mse[1]
  res[1,5]=mst[2];res[2,5]=msr[2];res[3,5]=mse[2]
  res[1,6]=f1[1];res[2,6]=f1[2]
  res[1,7]=f2[1];res[2,7]=f2[2]
  res[1,8]=qf(0.95,df2=(nm-n-m+1),df1=(n-1));res[2,8]=qf(0.95,df2=(nm-n-m+1),df1=(m-1))
  res=round(res,2)
  message(paste0(capture.output(res),collapse="\n"))
  message("--------------------------------------------------","\n")
  if(f1[1]>res[1,8])
  {
    message("we reject the null hypothesis","\n")
    message("so that there is a significant difference between different treatments","\n")
  }
  if(f2[1]<res[1,8])
  {
    message("we accept the null hypothesis","\n")
    message("so that there is no significant difference  between different treatments","\n")
  }
  if(res[1,8]>f1[1] && res[1,8]<f2[1])
  {
    message("we have no decision (we are indeterminant) regarding to the treatments","\n")
  }

  if(f1[2]>res[2,8])
  {
    message("we reject the null hypothesis","\n")
    message("so that there is a significant difference between different sectors","\n")
  }
  if(f2[2]<res[2,8])
  {
    message("we accept the null hypothesis","\n")
    message("so that there is no significant difference  between different sectors","\n")
  }
  if(res[2,8]>f1[2] && res[2,8]<f2[2])
  {
    message("we have no decision (we are indeterminant) regarding to the sectors","\n")
  }
}
