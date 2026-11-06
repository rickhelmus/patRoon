# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

wfHasSlot <- function(object, slotName, onlySetObj = FALSE)
{
    if (!is.null(slot(object, slotName)))
        return(TRUE)
    
    if (!object@dryRun || nrow(object@steps) == 0)
        return(FALSE)
    
    if (onlySetObj)
        return(slotName %chin% object@steps[!is.na(set)]$slotNameOut)
    return(slotName %chin% object@steps[is.na(set)]$slotNameOut)
}

wfHasScreeningFGroups <- function(obj, onlySetObj = FALSE)
{
    if (!wfHasSlot(obj, "fGroups", onlySetObj = onlySetObj))
        return(FALSE)
    
    if (obj@dryRun)
    {
        st <- obj@steps[func == "screenSuspectsP" & slotNameOut == "fGroups"]
        st <- if (onlySetObj)
            st[!is.na(set)]
        else
            st[is.na(set)]
        return(nrow(st) > 0)
    }
    
    return(isScreening(obj@fGroups))
}
