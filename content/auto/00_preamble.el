(TeX-add-style-hook
 "00_preamble"
 (lambda ()
   (TeX-run-style-hooks
    "pdfpages"
    "amsmath"
    "hyperref"
    "todonotes"))
 :latex)

