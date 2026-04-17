
ret = bettermc::mclapply(1:4, function(i) {
  if (i == 1L)
    stop(i)
  else if (i == 4L)
    system(paste0("kill ", Sys.getpid()))
  NULL
}, mc.allow.fatal = TRUE, mc.allow.error = TRUE, mc.preschedule = FALSE)


# retries
set.seed(456)
k=0
res <- bettermc::mclapply(1:2, function(i) {
	k=k+1
	if(k < 2) system(paste0("kill ", Sys.getpid()))
	k
  }, mc.retry = 3, mc.allow.fatal = TRUE, mc.cores = 5, mc.force.fork=T)




set.seed(456)
res <-
  bettermc::mclapply(1:20, function(i) {
    r <- runif(1)
    if (r < 0.25)
      system(paste0("kill ", Sys.getpid()))
    else if (r < 0.5)
      stop(i)
    else
      i
  }, mc.retry = 0, mc.cores = 10, mc.force.fork = TRUE, mc.allow.fatal = TRUE, mc.allow.error = TRUE)
