
######### CHANGE VALUE OF VARIABLE WITHIN A SLIM BLOCK  #########
udpate_block=function(block,...){
  args=match.call(expand.dots = TRUE)
  PARAMS_NAME=as.character(names(args[-(1:2)]))
  PARAMS_VALUE=as.character(sapply(args[-(1:2)],eval))
  for (i in 1: length(PARAMS_NAME)){
    block=gsub(PARAMS_NAME[i],PARAMS_VALUE[i],block)
  }
  block
 }

######### CONCATENATE SLIM BLOCKS #########
generateScript=function(...,quiet=FALSE){
  args=match.call(expand.dots = TRUE)
  args=args[as.character(names(args))!='quiet']
  BLOCKS=sapply(args[-(1)],function(x){as.character(eval(x),envir=parent.frame())})
  temp_script=sprintf('%s.slim',tempfile())
  if(!quiet){
    cat(paste0(BLOCKS,'\n',collapse='\n\n'))
  }
  paste0(BLOCKS,'\n',collapse='\n\n')
}

######### RUN SCRIPT AND READ OUTPUT #########
runScript=function(SCRIPT,
                    scriptFILE=sprintf('%s.slim',tempfile())
                    ){
  SCRIPT_clean=gsub('#','//',SCRIPT)
  SCRIPT_clean=gsub('\"\t\"','\"\\\\t\"',SCRIPT_clean)

  write(SCRIPT_clean,file=scriptFILE)
  SLIM_CMD=paste('slim',scriptFILE)
  out=system(SLIM_CMD,intern=TRUE)
  out
}
