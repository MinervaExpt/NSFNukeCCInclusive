SUBDIRS := src

all: $(SUBDIRS)

clean:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean; \
	done

arachne: src

$(SUBDIRS):
	$(MAKE) -C $@

.PHONY: all clean $(SUBDIRS)
