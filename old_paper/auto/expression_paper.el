(TeX-add-style-hook "expression_paper"
 (lambda ()
    (LaTeX-add-bibliographies
     "/home/ele/bibtex/master"
     "/home/ele/bibtex/master2")
    (LaTeX-add-environments
     "bmcformat")
    (TeX-add-symbols
     "includegraphic"
     "includegraphics")
    (TeX-run-style-hooks
     "longtable"
     "multirow"
     "inputenc"
     "utf8"
     "multicol"
     "ifthen"
     "url"
     "cite"
     "latex2e"
     "bmc_article10"
     "bmc_article"
     "10pt")))

