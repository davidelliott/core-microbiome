all.bac.rel.genus <- tax_glom(all.bac.rel,taxrank = "Genus")

kb.rel.genus <- tax_glom(kb.rel,taxrank = "Genus")
hb.rel.genus <- tax_glom(hb.rel,taxrank = "Genus")

all.bac.rel.genus.alt <- merge_phyloseq(kb.rel.genus,hb.rel.genus)
all.bac.rel.genus.alt.glom <- tax_glom(all.bac.rel.genus.alt,"Genus")

all.bac.rel.genus
all.bac.rel.genus.alt
all.bac.rel.genus.alt.glom

head(tax_table(kb.rel))
head(tax_table(hb.rel))




all.bac.rel.original <- all.bac.rel
all.bac.rel <- tax_glom(all.bac.rel.original, "Species")
