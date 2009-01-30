`generateLpatterns` <-
function()
#######################################################################
# all possible response patterns and/or difference patterns
#######################################################################
#
# CASE LIKERT
{

   # data: difference patterns
   diffsdat<-differences(get("dat",get("ENV",environment(patt.design))))

#print(head(diffsdat))

   # patterns
   ncatL<-get("ncatL",get("ENV",environment(patt.design)))
   patterns<-all_patterns(0,(ncatL-1))  # all possible likert patterns
   diffs<-differences(patterns)         # all possible diference patterns

   ## for the time being responses are centered around 0
   # shift data/patterns to 0:(ncatPC-1)
   #   diffs<-diffs-min(diffs)
   #   diffsdat<-diffsdat-min(diffsdat)

   if (get("blnReducecat",get("ENV",environment(patt.design)))) {
       # both next statements for -1/1 parameterisation (instead of 0/1)
       # small values in diff correspond to first obj preferred if blnRevert=FALSE
       diffs<-ifelse(diffs<0,1,ifelse(diffs>0,-1,0))
       diffsdat<-ifelse(diffsdat<0,1,ifelse(diffsdat>0,-1,0))
       assign("ncatL",2,envir=sys.frame(-1))
       assign("ncatPC",3,envir=get("ENV",environment(patt.design)))
   }

   # remove redundant patterns from diffs
   diffs<-unique(diffs, MARGIN=1)



   # convert diffs (patterns) to string
   dpattStr <- convert2strings(diffs)
   assign("datStr",convert2strings(diffsdat),envir=sys.frame(-1))


   # add to searchpath
   assign("dpattStr",dpattStr,envir=sys.frame(-1))       # character representation
   assign("npatt",length(dpattStr),envir=sys.frame(-1))  # number of unique possible patterns
   assign("diffs",diffs,envir=sys.frame(-1))             # numeric representation

   assign("blnUndec",TRUE,envir=sys.frame(-1))
}
