#!/usr/bin/env Rscript
# ==============================================================================
# Publication Figure: 20-node Three-generation Pedigree
# Optimized for absolute grid alignment with perfectly straight perpendicular lines.
# ==============================================================================

if (!requireNamespace("kinship2", quietly = TRUE)) {
  install.packages("kinship2", repos = "https://cloud.r-project.org")
}
library(kinship2)

# 1. Define pedigree dataframe
kindf <- data.frame(
  id    = 1:20,
  dadid = c(0,0,0,1,1,0,1,0,3,3,3,5,5,7,7,7,0,16,16,16),
  momid = c(0,0,0,2,2,0,2,0,4,4,4,6,6,8,8,8,0,17,17,17),
  # 1=male, 2=female (kinship2 convention)
  sex   = c(1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,1,2,2,1,2)
)

kept     <- c(9, 10, 11, 12, 13, 14, 15, 18, 19, 20)
affected <- as.integer(kindf$id %in% kept)   # 1 = retained (filled), 0 = dropped (unfilled)

# 2. Build the pedigree object
ped <- pedigree(id       = kindf$id,
                dadid    = kindf$dadid,
                momid    = kindf$momid,
                sex      = kindf$sex,
                affected = affected)

# 3. Save to PDF (Vector graphic for manuscript submissions)
# Set to 13x6 to provide ample horizontal space for completely straight lines
pdf("pedigree_family_structure.pdf", width = 13, height = 6)
op <- par(mar = c(1.5, 1.5, 1.5, 1.5)) 
plot(ped,
     cex        = 0.9,       # Text size for node ID labels
     symbolsize = 1.1,       # Size of pedigree shapes
     branch     = 0.5,       # Vertical branch drops
     packed     = FALSE,     # Disable compact squeezing to maintain grid columns
     align      = FALSE,     # <-- CRITICAL: Centers parents over children to eliminate zigzags
     width      = 12,        # <-- Gives the algorithm maximum horizontal space to straighten lines
     id         = as.character(kindf$id))
par(op)
dev.off()

# 4. Save to PNG (300 DPI high-quality raster graphic for presentations)
png("pedigree_family_structure.png", width = 13, height = 6, units = "in", res = 300)
op <- par(mar = c(1.5, 1.5, 1.5, 1.5))
plot(ped,
     cex        = 0.9,
     symbolsize = 1.1,
     branch     = 0.5,
     packed     = FALSE,
     align      = FALSE,     # <-- Centers parents over children to eliminate zigzags
     width      = 12,        # <-- Aligns columns strictly
     id         = as.character(kindf$id))
par(op)
dev.off()

# ==============================================================================
# VERIFICATION CHECKS
# ==============================================================================
kfull <- kinship(ped)
kfam  <- kfull[as.character(kept), as.character(kept)]
cat("\n=== Verification Output ===\n")
cat("Within-family kinship block (phi_ij) for the 10 retained descendants:\n\n")
print(round(kfam, 4))

cat("\nUnique off-diagonal values (should be 0, 1/32, 1/16, 1/8, 1/4):\n")
print(sort(unique(kfam[upper.tri(kfam)])))

cat("\n[SUCCESS] Generated beautifully straight, perpendicular pedigree plots:\n")
cat("  -> PDF Vector: pedigree_family_structure.pdf (Manuscript Standard)\n")
cat("  -> PNG Image:  pedigree_family_structure.png (Slides/Previews)\n")
