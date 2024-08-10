# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = python -msphinx
SPHINXPROJ    = BMO
SOURCEDIR     = .
BUILDDIR      = _build
SKIP_BUILD_DOC =  # Set this to '1' to skip `python ../build_doc.py`

# Put it first so that "make" without argument is like "make help".
help:
	#@python ../build_doc.py
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile clean all local clean_tmp


# Collect the targets as tags
TAGS = $(subst html,html,$(subst singlehtml,singlehtml,$(subst latexpdf,pdf,$(subst epub,epub,$(subst man,man,$(filter-out Makefile,$(MAKECMDGOALS)))))))


# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@if [ -z "$(SKIP_BUILD_DOC)" ]; then python ../build_doc.py; fi
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O) $(addprefix -t ,$(TAGS))

clean:
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O) $(addprefix -t ,$(TAGS))
	@find "../../python/" -type d -name "__pycache__" -print0 2>/dev/null | xargs -0 rm -rf

all:
	@if [ -z "$(SKIP_BUILD_DOC)" ]; then python ../build_doc.py; fi
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O) -t html
	@$(SPHINXBUILD) -M singlehtml "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O) -t singlehtml
	@$(SPHINXBUILD) -M latexpdf "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O) -t pdf
	@$(MAKE) clean_tmp

local:
	@if [ -z "$(SKIP_BUILD_DOC)" ]; then python ../build_doc.py; fi
	@$(SPHINXBUILD) -M singlehtml "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O) -t singlehtml
	@$(SPHINXBUILD) -M latexpdf "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O) -t pdf
	@$(MAKE) clean_tmp

clean_tmp:
	@find "$(BUILDDIR)"/epub/ -type f ! -name "*.epub" -print0 2>/dev/null | xargs -0 rm -f
	@find "$(BUILDDIR)"/epub/ -type d -empty -print0 2>/dev/null | xargs -0 rmdir
	@find "$(BUILDDIR)"/latex/ -type f ! -name "*.pdf" -print0 2>/dev/null | xargs -0 rm -f
	@find "$(BUILDDIR)"/latex/ -type d -empty -print0 2>/dev/null | xargs -0 rmdir
	@rm -rf "$(BUILDDIR)"/doctrees
	@find "../../python/" -type d -name "__pycache__" -print0 2>/dev/null | xargs -0 rm -fr
