`generateRpatterns` <-
function()
#######################################################################
# all possible response patterns and/or difference patterns
#######################################################################
#
# CASE RANKINGS
{

   # data: difference patterns
   diffsdat<-differences(get("dat",get("ENV",environment(patt.design))))

   # patterns
   nobj<-get("nobj",get("ENV",environment(patt.design)))
   patterns<-permutations(nobj)         # all possible ranks patterns
   diffs<-differences(patterns)         # all possible diference patterns

   # for rankings reduced categories
   # -> same result as ordinal categories
   # -> alway reduce categories

   # both next statements for -1/1 parameterisation (instead of 0/1)
   # small values in diff correspond to first obj preferred if blnRevert=FALSE
   diffs<-ifelse(diffs<0,1,-1)   # rankings never undecided
   diffsdat<-ifelse(diffsdat<0,1,-1)
   assign("ncatPC",2,envir=sys.frame(-1))
   assign("blnUndec",FALSE,envir=sys.frame(-1))

   # convert diffs (patterns) to string
   dpattStr <- convert2strings(diffs)
   assign("datStr",convert2strings(diffsdat),envir=sys.frame(-1))


   # add to searchpath
   assign("dpattStr",dpattStr,envir=sys.frame(-1))       # character representation
   assign("npatt",length(dpattStr),envir=sys.frame(-1))  # number of unique possible patterns
   assign("diffs",diffs,envir=sys.frame(-1))             # numeric representation
}
