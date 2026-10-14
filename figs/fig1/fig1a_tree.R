library(ape)

# read the IQ-TREE species tree and root it on the mayflies
tr  <- read.tree("../../palaeoptera_synteny/pairwise_synteny/tree/iqtree/tree.treefile")
eph <- c("Baetis_atlanticus","Baetis_rhodani","Cloeon_dipterum","Procloeon_bifidum",
         "Ecdyonurus_torrentis","Parafronurus_youi","Ephemera_orientalis",
         "Ephemera_vulgata","Serratella_ignita","Siphlonurus_alternatus")
tr  <- ladderize(root(tr, outgroup = eph, resolve.root = TRUE))

# put the root halfway between the orders, drop the "Root" node label, tidy tip names
re  <- which(tr$edge[, 1] == Ntip(tr) + 1)          # the two branches at the root
tr$edge.length[re] <- sum(tr$edge.length[re]) / 2
tr$node.label[tr$node.label == "Root"] <- ""
tr$tip.label <- gsub("_", " ", tr$tip.label)

# plot with aligned tip labels and a scale bar
svg("palaeoptera.tree.svg", width = 12, height = 33.23)
plot(tr, font = 3, cex = 3, label.offset = 0.01, no.margin = TRUE,
     align.tip.label = TRUE,
     edge.width = 3,
     x.lim = c(0, 1.5 * max(node.depth.edgelength(tr))))   # room for the aligned labels
add.scale.bar(length = 0.2, lwd = 3)              # substitutions/site
dev.off()

# print the tip order, top to bottom as plotted
is_tip <- tr$edge[,2] <= length(tr$tip.label)
ordered_tips <- tr$edge[is_tip, 2]
rev(tr$tip.label[ordered_tips])
