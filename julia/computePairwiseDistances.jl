#computes pairwise distances between sequences in a fasta file

# Reading a fasta file
function readFasta(file)
  seqs = AbstractString[]
  nams = AbstractString[]
  seq = ""
  seqname = ""
  open(file) do f
    for i in enumerate(eachline(f))
      if startswith(i[2], '>')
        if seq != "" && seqname != ""
          push!(seqs,seq)
          push!(nams,seqname)
        end
        seq=""
        liste=split(i[2], ":")
        freq = strip(liste[2])
        if ( freq != "0")
          seqname = replace(liste[1], ">","")
        else
          seqname = ""
        end
      else
        seq = seq * strip(i[2])
      end
    end
  end
  # last sequence
  if seq != "" && seqname != ""
      push!(seqs, seq)
      push!(nams, seqname)
  end
  println("Number of non-0 haplotypes: ", length(seqs) )
  return (nams, seqs)
end

# Computing the pairwise differences between two sequences
function computePairwiseDifferences( x, y )
  dif = 0
  for k = 1:length(x)
    xk = x[k]
    yk = y[k]
    if xk != yk && xk != 'N' && yk != 'N'
      dif += 1
    end
  end
  return dif
end


function computeAndOutputPairwiseDifferences(out, nams, seqs)
  open(out, "w") do f
    write(f, "from\tto\tdist\n")
    for seqi = 1:(length(seqs)-1)
      namei = nams[seqi]
      for seqj = (seqi+1):length(seqs)
        dif = computePairwiseDifferences( seqs[seqi], seqs[seqj] )
        write(f, namei*"\t"*nams[seqj]*"\t"*string(dif), "\n")
      end
    end
  end
#    if ( (seqi % 10) == 0):
#        print ("seqi: "+str(seqi))
end


file=ARGS[1]
out=ARGS[2]

function main()
  ( nams, seqs ) = readFasta(file)
  computeAndOutputPairwiseDifferences(out, nams, seqs)
end

# The following needs to be commented out for compilation into a binary
main()
