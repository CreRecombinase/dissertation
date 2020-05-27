((org-mode . ((org-publish-project-alist .
                                         (
                                          ("dissertation-org" . (:base-directory "/home/nwknoblauch/Dropbox/Repos/dissertation/org"
                                                                                 :publishing-directory "/home/nwknoblauch/Dropbox/Repos/dissertation/docs"
                                                                                 :publishing-function org-html-publish-to-html
                                                                                 :makeindex t
                                                                                 :auto-sitemap t
                                                                                 :sitemap-filename "index.org"
                                                                                 :sitemap-title "Home"
                                                                                 :sitemap-sort-files anti-chronologically
                                                                                 :sitemap-file-entry-format "%d - %t"
                                                                                 :sitemap-function org-publish-org-sitemap
                                                                                 :exclude "setup.org"
                                                                                 ))
                                          ("dissertation-fig" . (:base-directory "/home/nwknoblauch/Dropbox/Repos/dissertation/org"
                                                                                 :publishing-directory "/home/nwknoblauch/Dropbox/Repos/dissertation/docs"
                                                                                 :base-extension "png\\|\\|svg"
                                                                                 :publishing-function org-publish-attachment
                                                                                 :recursive t
                                                                                 ))
                                          ("dissertation" . (:components ("dissertation-org" "dissertation-fig") )))
                                         )
              (reftex-default-bibliography '("~/Dropbox/Repos/dissertation/bibliography/references.bib"))
              (org-ref-bibliography-notes "~/Dropbox/Repos/dissertation/bibliography/notes.org")
              (org-ref-pdf-directory "~/Dropbox/Repos/dissertation/bibliography/bibtex-pdfs/")

              )))
