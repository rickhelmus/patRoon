getRDepsTab <- function(format)
{
    rdpath <- tempfile(fileext = ".R")
    download.file("https://rickhelmus.github.io/patRoonDeps/utils/Rdeps.R", rdpath)
    rdenv <- new.env()
    source(rdpath, local = rdenv)
    
    depsList <- rdenv$getRDependencies("master", os = NULL, onlyPDeps = FALSE, withInternal = FALSE, flatten = TRUE)
    
    availPD <- available.packages(repos = "https://rickhelmus.github.io/patRoonDeps", type = "binary")
    availRU <- available.packages(repos = "https://rickhelmus.r-universe.dev")
    
    ghRepos <- function(name, dep) paste0(dep$user, "/", if (!is.null(dep[["repos"]])) dep$repos else name)
    ghDepURL <- function(name, dep)
    {
        repos <- ghRepos(name, dep)
        if (!is.null(dep[["pkgroot"]]))
            repos <- paste0(repos, "/", dep$pkgroot)
        return(sprintf("https://github.com/%s", repos))
    }
    depComments <- function(dep)
    {
        # parentDep, OS, mandatory, tag, branch, patRoonDeps
        ret <- character()
        add <- function(val, what) if (!is.null(dep[[val]])) ret <<- c(ret, what)
        add("parentDep", paste("Dependency of", dep$parentDep))
        add("os", paste("Only for", dep$os))
        add("mandatory", "Mandatory")
        return(paste0(ret, collapse = "<br>"))
    }
    depInstPD <- function(name, dep)
    {
        if (!name %in% rownames(availPD))
            return("N/A")
        return(sprintf("install.packages('%s', repos = 'https://rickhelmus.github.io/patRoonDeps', type = 'binary')", name))
    }
    depInstRU <- function(name, dep)
    {
        if (!name %in% rownames(availRU))
            return("N/A")
        return(sprintf("install.packages('%s', repos = 'https://rickhelmus.r-universe.dev')", name))
        
    }
    depInstReg <- function(name, dep)
    {
        if (dep$type == "cran")
            return(sprintf("`install.packages('%s')`", name))
        if (dep$type == "bioc")
            return(sprintf("`BiocManager::install('%s')`", name))
        
        # else GitHub
        
        repos <- ghRepos(name, dep)
        for (ref in c("branch", "tag", "commit"))
        {
            if (!is.null(dep[[ref]]))
            {
                repos <- paste0(repos, "@", dep[[ref]])
                break
            }
        }
        return(sprintf("`remotes::install_github('%s')`", repos))
    }
    
    tab <- data.table::data.table(name = names(depsList), type = lapply(depsList, "[[", "type"))
    tab[, url := data.table::fcase(type == "cran", sprintf("https://cran.r-project.org/web/packages/%s/index.html", name),
                                   type == "bioc", sprintf("https://bioconductor.org/packages/release/bioc/html/%s.html", name),
                                   type == "gh", mapply(name, depsList[name], FUN = ghDepURL))]
    if (format == "html")
        tab[, package := sprintf("<a href='%s'>%s</a>", url, name)]
    else # latex
        tab[, package := sprintf("\\href{%s}{%s}", url, name)]
    tab[, comments := sapply(depsList, depComments)]
    # tab[, installPD := mapply(name, depsList[name], FUN = depInstPD)]
    # tab[, installRU := mapply(name, depsList[name], FUN = depInstRU)]
    tab[, patRoonDeps := data.table::fifelse(name %in% rownames(availPD), "yes", "no")]
    tab[, `r-universe` := data.table::fifelse(name %in% rownames(availRU), "yes", "no")]
    tab[, `regular installation` := mapply(name, depsList[name], FUN = depInstReg)]
    
    tab[, c("name", "type", "url") := NULL]
    
    return(tab)
}

