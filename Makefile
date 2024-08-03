
out/%.html: src/%.md
	pandoc --standalone --toc --output $@ $<
