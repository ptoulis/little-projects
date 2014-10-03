output.files = list.files("odyssey", full.names=T, pattern="rda")
A <- c()
for(file in output.files) {
  load(file)
  A = c(A, out.coverage)
}

print("Coverage statistics")
print(summary(A))

