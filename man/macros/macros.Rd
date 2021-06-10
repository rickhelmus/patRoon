\newcommand{\addCitations}{\Sexpr[results=text,echo=FALSE,strip.white=FALSE]{format(citation("#1")[[#2]], "textVersion")}}

\newcommand{\setsWF}{(\strong{sets workflow})}
\newcommand{\setsWFClass}{The \code{#1} class is applicable for \link[=sets-workflow]{sets workflows}. This class is derived from \code{#2} and therefore largely follows the same user interface.}
\newcommand{\setsWFUnset}{\item \code{unset} Converts the object data for a specified set into a 'non-set' object (\code{#1}), which allows it to be used in 'regular' workflows. #2}
\newcommand{\setsWFNewMethodsFeat}{The following methods are specifically defined for sets workflows: \itemize{ \item \code{sets} Returns the set names for this object. \setsWFUnset{#1}{#2} }}
\newcommand{\setsWFNewMethodsSO}{The following methods are specifically defined for sets workflows: \itemize{\item All the methods from base class \code{\link{workflowStepSet}}. \setsWFUnset{#1}{#2} }}
\newcommand{\setsWFNewMethodsSOExtra}{The following methods are specifically defined for sets workflows: \itemize{\item All the methods from base class \code{\link{workflowStepSet}}. \setsWFUnset{#1}{#2} #3 }}
\newcommand{\setsWFChangedMethods}{The following methods are changed or with new functionality: \itemize{#1}}

