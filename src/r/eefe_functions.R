# Color mutations
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  Non_Synonymous = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#A6CEE3", col = NA))
  },
  Small_indel = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#FB9A99", col = NA))
  },
  Intergenic_snp = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = '#FF7F00', col = NA))
  },
 Stop = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = '#6A3D9A', col = NA))
  },
  Synonymous = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#FDBF6F", col = NA))
  },
  Large_deletion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = '#E31A1C', col = NA))
  },
  Large_amplification = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#1F78B4", col = NA))
  }
)
col = c("Non_Synonymous" = "#A6CEE3", "Small_indel" = "#FB9A99", "Intergenic_snp" = '#FF7F00',
        'Stop'='#6A3D9A','Synonymous'='#FDBF6F','Large_deletion'='#E31A1C',"Large_amplification"="#1F78B4")