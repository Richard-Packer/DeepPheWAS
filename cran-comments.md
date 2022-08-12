## Test Enviroments
* Local CentOS Linux 7 (core), R 4.2.1
* Local Windows 10, R 4.2.1
* Local OS X, R 4.2.1
* Ubuntu Linux 20.04.1 LTS (on R-hub), R-release
* Fedora Linux (on R-hub), R-devel
* Windows Server 2022 (on R-hub), R-devel 

## R CMD check results
There were no ERRORs or WARNINGs. 

There were two NOTES found on Windows server 2022. The first note is also found on Ubuntu Linux 
and Fedora Linux.

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Richard packer <richard.packer@leicester.ac.uk>'

New submission

Possibly misspelled words in DESCRIPTION:
  Phenome (3:13)
  Phenotypes (2:59)
  phenome (8:54)
  phenotypes (8:31)


Possibly misspelled words in DESCRIPTION:
    Phenome (3:13)
    Phenotypes (2:59)
    phenome (8:54)
    phenotypes (8:31)

* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'
```
The first note none of these words are misspelt.
The latter NOTE As noted in [R-hub issue #503](https://github.com/r-hub/rhub/issues/503), this could be due to a bug/crash in MiKTeX and can likely be ignored.

## Downstream dependencies
There are currently no down stream dependencies

---

This is the first version being uploaded.


Thanks,
Richard Packer
