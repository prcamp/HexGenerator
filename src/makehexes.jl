Scl = 1
const α = Scl*(sqrt(3)/2)
const β = Scl*(1/2)


c(i,N) = (i>N)? 2*N-i : i

# compute the number of vertices in a given row
nv(i,numsides) = i+numsides
numvertsthisrow(i,numsides) = (i <= numsides + 1)? nv(i,numsides) : nv(c(i,numsides+1),numsides)

# compute the number of faces in a given row
nf(i,numsides) = 2*numsides + (2*i-1)
numfacesthisrow(i,numsides) = (i <= numsides)? nf(i,numsides) : nf(c(i-1,numsides),numsides)


# Given coordinates for a face, get coordinates for the three vertices:
function face2verts(i,j,numsides)
  verts = Vector{Tuple{Int,Int}}()
  k = (iseven(j))? round(Int,j/2) : round(Int,(j+1)/2)
  if i <= numsides
    if !iseven(j)
      verts = [(i,k),(i+1,k),(i+1,k+1)]
    else
      verts = [(i,k),(i+1,k+1),(i,k+1)]
    end
  else
    if !iseven(j)
      verts = [(i,k),(i+1,k),(i,k+1)]
    else
      verts = [(i,k+1),(i+1,k),(i+1,k+1)]
    end
  end
  return verts
end

jhat(i,j,numsides)  = (i <= numsides)? j : j+2*(i-numsides)-1
ijhat(i,j,numsides) = (i <= numsides)? j : j+1-2*(i-numsides)

function faceneighbors(i,j,numsides)
  faceneighbors = Vector{Tuple{Int,Int}}()
  jh = jhat(i,j,numsides)
  nftr = numfacesthisrow(i,numsides)
  (j > 1) && push!(faceneighbors,(i,j-1))
  (j < nftr) && push!(faceneighbors,(i,j+1))
  if !iseven(jh) && (i < 2*numsides)
    push!(faceneighbors,(i+1,ijhat(i+1,jh+1,numsides)))
  elseif (i > 1)
    push!(faceneighbors,(i-1,ijhat(i-1,jh-1,numsides)))
  end
end

isinterior(i,j,numsides) = !(i == 1 || i == 2*numsides+1 || j==1 || j==numvertsthisrow(i,numsides))

function vertneighbors(i,j,numsides)
  nbors = Vector{Tuple{Int,Int}}()
  if isinterior(i,j,numsides)
    if i <= numsides+1
      push!(nbors,(i-1,j-1))
      push!(nbors,(i-1,j))
    else
      push!(nbors,(i-1,j))
      push!(nbors,(i-1,j+1))
    end
    push!(nbors,(i,j-1))
    push!(nbors,(i,j+1))
    if i < numsides+1
      push!(nbors,(i+1,j))
      push!(nbors,(i+1,j+1))
    else
      push!(nbors,(i+1,j-1))
      push!(nbors,(i+1,j))
    end
  end
  return nbors
end

function computediffs(vertvalues)
  diffs = Dict{Tuple{Int,Int},Int}()
  numsides = length(vertvalues[1])-1
  for (i,vv) in enumerate(vertvalues)
    for (j,v) in enumerate(vv)
      vertnbrs = vertneighbors(i,j,numsides)
      if !isempty(vertnbrs)
        nbrsum = 0
        for ind in vertnbrs
          nbrsum += vertvalues[ind[1]][ind[2]]
        end
        nbrsum -= 6*vertvalues[i][j]
        diffs[(i,j)]=nbrsum
      end
    end
  end
  diffs
end

function compute_diffs_and_flippable(vertvalues)
  diffs = computediffs(vertvalues)
  flippable = Set{Tuple{Int64,Int64}}()
  negflippable = Set{Tuple{Int64,Int64}}()
  for (k,v) in diffs
    if (v==9) || (v==-9)
      push!(flippable,k)
    end
  end
  return diffs,flippable,negflippable
end

isflippable(x) = (x==9) || (x==-9)
isnegflippable(x) = (x==-9)

function ind2coords(i,j,numsides,scl=1)
  #α = (sqrt(3)/2)
  #β = (1/2)
  # XXXXXXXXXXXX
  initialpoint = [-scl*α,scl*β] # [-3*scl*α,scl*β] #
  ur = [scl*α,scl*β]/numsides
  d = [0,-2*scl*β]/numsides
  dr = [scl*α,-scl*β]/numsides
  if i <= numsides+1
    return  initialpoint + (i-1)*d + (j-1)*ur
  else
    return initialpoint + numsides*d + (i-numsides-1)*dr + (j-1)*ur
  end
end

function faceorientation(vertvalues,i,j,numsides)
  vals = Vector{Int}()
  for v in face2verts(i,j,numsides)
    push!(vals,vertvalues[v[1]][v[2]])
  end
  return sortperm(vals)
end

function makeface(vals,vertcoords,numsides,macroscl=1)
  scl = macroscl/numsides
  const type1 = [2,1,3]
  const type2 = [1,3,2]
  const type3 = [3,2,1]
  fo = sortperm(vals)
  ft = 0
  if fo == type1
    ft = 1
    center = .5*(vertcoords[2,:]+vertcoords[3,:])
  elseif fo == type2
    ft = 2
    center = .5*(vertcoords[1,:]+vertcoords[2,:])
  elseif fo == type3
    ft = 3
    center = .5*(vertcoords[1,:]+vertcoords[3,:])
  else
    @show vals
    Base.error("The orientation of that face is not one of the types")
  end
  fc = facecoords(ft, center, scl)
  return fc,ft
end

function center(face)
  1/2*(face[1]+face[3])
end


function makefaces(vertvalues,numsides,scl=1,center=(0,0))
  numfaces = numsides^2*3
  facescoords = Curves(numfaces)#Vector{Array{Float64}}()#
  facetypes = zeros(Int,numfaces)
  cent = [center[1] center[2]]
  #    cents = Vector{Tuple{Float64,Float64}}()
  facecnt = 0
  vert2facedict = Dict{Tuple{Int,Int},Set{Int}}()
  for i in 1:2*numsides
    nftr = numfacesthisrow(i,numsides)
    for j in 1:nftr
      if ((i <= numsides) && !iseven(j)) || ((i > numsides) && iseven(j))
        faceverts = face2verts(i, j, numsides)
        vals = Vector{Int}(3)
        vertcoords = Array{Float64}(3,2)
        facecnt += 1
        for k in 1:3
          v = faceverts[k]
          (!haskey(vert2facedict,v)) && (vert2facedict[v] = Set{Int}())
          push!(vert2facedict[v],facecnt)
          vertcoords[k,:] = (ind2coords(v[1], v[2], numsides,scl)[:]')#+cent
          vals[k] = vertvalues[v[1]][v[2]]
        end
        fc,ft = makeface(vals,vertcoords,numsides,scl)
        face = Curve(5)
        for k in 1:5
          face[k] = Point(fc[k,1],fc[k,2])
        end
        facescoords[facecnt] = face
        facetypes[facecnt] = ft
        #                push!(cents,(cntr[1],cntr[2]))
        #face = mapslices(f->Point(f[1],f[2]),fc,2)
        #push!(facescoords,fc)
        #push!(facetypes,ft)
      end
    end
  end
  facecentx = .5*(maxx(facescoords)+minx(facescoords))
  facecenty = .5*(maxy(facescoords)+miny(facescoords))
  facecent = Point(facecentx,facecenty)
  cent = Point(cent[1],cent[2])
  for (i,f) in enumerate(facescoords)
    f += (-1)*facecent
    f += cent
    facescoords[i] = f
  end
  return facescoords,facetypes,vert2facedict
end

function facecoords(facetype,center=[0,0],mcrscl=1)
  fα = mcrscl*α
  fβ = mcrscl*β
  if facetype == 1
    return 1/2*[-fα -fβ; -fα 3*fβ; fα fβ; fα -3*fβ;-fα -fβ]+repmat(center[:]',5,1)
  elseif facetype == 2
    return [ -fα 0 ; 0 fβ ; fα 0 ; 0 -fβ; -fα 0]+repmat(center[:]' ,5,1)
  elseif facetype == 3
    return 1/2*[fα -fβ; fα 3*fβ ; -fα fβ ; -fα -3*fβ; fα -fβ]+repmat(center[:]' ,5,1)
  end
end

function geninitial0(numsides)
  numrows = 2*numsides+1
  nsp1 = numsides+1
  vertvalues = Vector{Vector{Int}}()
  for i in 1:numrows
    nvtr = numvertsthisrow(i,numsides) #num verts this row
    if i <= numsides+1
      valsthisrow1 = (2*nsp1)-i-collect(1:nsp1)
      valsthisrow2 = valsthisrow1[end]+2*collect(1:(nvtr-nsp1))
      valsthisrow = vcat(valsthisrow1,valsthisrow2)
      push!(vertvalues, valsthisrow)
    else
      icomp = 2*nsp1 - i
      valsthisrow1 = (2*nsp1)-icomp-collect(1:nsp1)
      valsthisrow2 = valsthisrow1[end]+2*collect(1:(nvtr-nsp1))
      valsthisrow = vcat(valsthisrow1,valsthisrow2)
      push!(vertvalues,valsthisrow)
    end
  end
  diffs = Dict{Tuple{Int,Int},Int}()
  flippable = Set{Tuple{Int,Int}}()
  negflippable = Set{Tuple{Int,Int}}()
  push!(flippable,(nsp1,nsp1))
  for i in 2:numrows-1
    for j in 2:numvertsthisrow(i,numsides)-1
      diffs[(i,j)]=0
      if (i == nsp1 && j < nsp1) || (i != nsp1 && j == nsp1)
        diffs[(i,j)]=6
      elseif i==nsp1 && j == nsp1
        diffs[(i,j)]=9
      end
    end
  end
  return vertvalues, diffs, flippable,negflippable
end

jverthat(i,j,numsides) = (i <= numsides+1)? j : j + (i-numsides-1)
ijverthat(i,j,numsides) = (i <= numsides+1)? j : j-(i-numsides-1)

function flipvert(fv::Tuple{Int,Int},vertvalues, diffs, flippable,negflippable,numsides)
  dir = diffs[fv]
  val = vertvalues[fv[1]][fv[2]]
  d = (sign(dir)>0)? 3 : -3
  vertvalues[fv[1]][fv[2]] += d
  diffs[fv] -= 6*d
  (!isflippable(diffs[fv])) && delete!(flippable,fv)
  (!isnegflippable(diffs[fv])) && delete!(negflippable,fv)
  for n in vertneighbors(fv[1],fv[2],numsides)
    if isinterior(n[1], n[2], numsides)
      diffs[n] += d
      if isflippable(diffs[n])
        push!(flippable,n)
      elseif isnegflippable(diffs[n])
        push!(negflippable,n)
      else
        delete!(flippable,n)
        delete!(negflippable,n)
      end
    end
  end
  return vertvalues, diffs, flippable,negflippable
end

function geninitial(numsides)
  nsp1 = numsides+1
  vertvalues, diffs, flippable,negflippable = geninitial0(numsides)
  # these are in jverthat coords
  flipped = Set{Tuple{Int,Int}}()
  toflip = Set{Tuple{Int,Int}}([(nsp1,nsp1)])
  stopif(F) = (2,nsp1) ∈ F
  while !stopif(flipped)
    tf = map(v->v,toflip)
    for vert in tf
      i = vert[1]
      jh = vert[2]
      j = ijverthat(i,jh,numsides)
      fv = (i,j)
      vertvalues,diffs,flippable,negflippable = flipvert(fv,vertvalues,diffs,flippable,negflippable,numsides)
      push!(toflip,(i-1,jh))
      push!(toflip,(i,jh-1))
      push!(toflip,(i+1,jh+1))
      delete!(toflip,(i,jh))
      push!(flipped,(i,jh))
    end
  end
  return vertvalues, diffs, flippable,negflippable
end

function geninitial_alternate(numsides)
  allverts,~,~,~ = geninitial0(1)
  for i = 2:numsides
    (iseven(i)) && (map!(v->v+3,allverts))
    allVerts = Vector{Vector{Int}}()
    vv,~,~,~ = geninitial0(i)
    for i in 1:length(allverts)
      allverts[i] = vcat(vv[i+1][1],allverts[i],vv[i+1][end])
    end
    start_vec = Vector{Vector{Int}}()
    push!(start_vec,vv[1])
    end_vec = Vector{Vector{Int}}()
    push!(end_vec,vv[end])
    allverts = vcat(start_vec,allverts,end_vec)
  end
  dfs,flp,nflp = compute_diffs_and_flippable(allverts)
  return allverts,dfs,flp,nflp
end

function randflip(vertvalues,diffs,flippable,negflippable,numsides)
  (!isempty(flippable))? fv = rand(map(i->i,flippable)) : fv = rand(map(i->i,negflippable))
  vertvalues,diffs,flippable,negflippable = flipvert(fv,vertvalues,diffs,flippable,negflippable,numsides)
  return vertvalues,diffs,flippable,negflippable
end

function genhex(numsides,numflips,scl=1,center=(0,0))
  # initial = "down"
  vertvalues, diffs, flippable, negflippable = geninitial0(numsides)
  # (initial=="down") && ()
  # (initial=="mid") && (vertvalues, diffs, flippable, negflippable = geninitial(numsides))
  for i in 1:numflips
    vertvalues,diffs,flippable,negflippable = randflip(vertvalues,diffs,flippable,negflippable,numsides)
  end
  faces,facetypes,vert2facedict = makefaces(vertvalues,numsides,scl,center)
  return faces,facetypes,vert2facedict,vertvalues,diffs,flippable,negflippable
end



# function genhex(numsides,numflips,scl=1,center=(0,0))
#   vertvalues, diffs, flippable, negflippable = geninitial(numsides)
#   for i in 1:numflips
#     vertvalues,diffs,flippable,negflippable = randflip(vertvalues,diffs,flippable,negflippable,numsides)
#   end
#   faces,facetypes = makefaces(vertvalues,numsides,scl,center)
#   return faces,facetypes,vertvalues,diffs,flippable,negflippable
# end
#
# function genrecurhex(numsides::Vector{Int},numflips,scl,depth=2,center=(0,0),allfaces = Curves(),depths=Int[],facetypes = Int[];probrecur=.5)
#   (depth==0) && (return allfaces,depths,facetypes)
#   vertvalues, diffs, flippable, negflippable = geninitial_alternate(numsides)
#   for i in 1:numflips
#     vertvalues,diffs,flippable,negflippable = randflip(vertvalues,diffs,flippable,negflippable,numsides)
#   end
#
#   faces0,fts0,v2f = makefaces(vertvalues,numsides,scl,center)
#   # @show faces
#   # println("\n\n\n")
#   # @show fts
#   # println("\n\n\n")
#   # @show v2f
#   centers = Vector{Tuple{Float64,Float64}}()
#   scales = Float64[]
#
#   # Now go through and make sure they don't overlap:
#   toflip = Set{Tuple{Int,Int}}()
#   faces2delete = Set(Int[])
#   if depth > 1
#     for f in flippable
#       if rand() < probrecur
#         # check to see if it overlaps with a vertex in toflip:
#         addtoflip = true
#         for tf in toflip
#           if !isempty(intersect(v2f[f],v2f[tf]))
#             addtoflip = false
#           end
#         end
#         @show addtoflip
#         if addtoflip
#           push!(toflip,f)
#           for f2d in v2f[f]
#             push!(faces2delete,f2d)
#           end
#           #vcat(faces2delete,map(identity,v2f[f]))
#           @show (v2f[f],toflip,faces2delete)
#           #centvec = ind2coords(f[1],f[2],numsides,scl)#(0,0)#compute center
#           boxtodelete = Curves()
#           for fcidx in v2f[f]
#             push!(boxtodelete,faces0[fcidx])
#           end
#           centvecx = .5*(maxx(boxtodelete)+minx(boxtodelete))
#           centvecy = .5*(maxy(boxtodelete)+miny(boxtodelete))
#           cent = (centvecx,centvecy)
#           push!(centers,cent)
#           scale = (1/numsides)*scl
#           push!(scales,scale)
#         end
#       end
#     end
#
#     # Now cull the appropriate faces from faces0:
#     faces = Curves()
#     fts = Int[]#fts0#Int[]
#     # @show faces2delete
#     for (i,f) in enumerate(faces0)
#       (!in(i,faces2delete)) && (push!(faces,f);push!(fts,fts0[i]))
#     end
#   else
#     faces = faces0
#     fts = fts0
#   end
#   allfaces=vcat(allfaces,faces)
#   facetypes = vcat(facetypes,fts)
#   depths = vcat(depths,depth*(ones(Int,length(faces))))
#   # @show numsides
#   # @show numflips
#   # @show scale
#   # @show centers
#   # @show depth
#   # @show allfaces
#   # @show depths
#   # @show facetypes
#   for (ctr,scl) in zip(centers,scales)
#     a,d,ft = genrecurhex(numsides,numflips,scl,depth-1,ctr,allfaces,depths,facetypes)
#     allfaces = vcat(allfaces,a)
#     depths = vcat(depths,d)#epth*(ones(Int,length(a))))
#     facetypes = vcat(facetypes,ft)
#   end
#   # @show allfaces
#   # @show depths
#   # @show facetypes
#   #FC,DP,FT =
#   #else
#   #    cent = (0,0)
#   #       genrecurhex(numsides,numsides,(1/numsides)*scl,cent,0,allfaces,depths,facetypes)
#   return allfaces,depths,facetypes
# end
#

function genrecurhex(numsides,numflips,scl,depth=2,center=(0,0),allfaces = Curves(),depths=Int[],facetypes = Int[]; probrecur=.5, initial_config = 0)

  (depth==0) && (return allfaces,depths,facetypes)
  G = geninitial
  if initial_config == 0
    G = geninitial
  elseif initial_config == 1
    G = geninitial0
  else initial_config == 2
    G = geninitial_alternate
  end
  vertvalues, diffs, flippable, negflippable = G(numsides)
  for i in 1:numflips
    vertvalues,diffs,flippable,negflippable = randflip(vertvalues,diffs,flippable,negflippable,numsides)
  end

  faces0,fts0,v2f = makefaces(vertvalues,numsides,scl,center)
  # @show faces
  # println("\n\n\n")
  # @show fts
  # println("\n\n\n")
  # @show v2f
  centers = Vector{Tuple{Float64,Float64}}()
  scales = Float64[]

  # Now go through and make sure they don't overlap:
  toflip = Set{Tuple{Int,Int}}()
  faces2delete = Set(Int[])
  if depth > 1
    for f in flippable
      if rand() < probrecur
        # check to see if it overlaps with a vertex in toflip:
        addtoflip = true
        for tf in toflip
          if !isempty(intersect(v2f[f],v2f[tf]))
            addtoflip = false
          end
        end
        @show addtoflip
        if addtoflip
          push!(toflip,f)
          for f2d in v2f[f]
            push!(faces2delete,f2d)
          end
          #vcat(faces2delete,map(identity,v2f[f]))
          @show (v2f[f],toflip,faces2delete)
          #centvec = ind2coords(f[1],f[2],numsides,scl)#(0,0)#compute center
          boxtodelete = Curves()
          for fcidx in v2f[f]
            push!(boxtodelete,faces0[fcidx])
          end
          centvecx = .5*(maxx(boxtodelete)+minx(boxtodelete))
          centvecy = .5*(maxy(boxtodelete)+miny(boxtodelete))
          cent = (centvecx,centvecy)
          push!(centers,cent)
          scale = (1/numsides)*scl
          push!(scales,scale)
        end
      end
    end

    # Now cull the appropriate faces from faces0:
    faces = Curves()
    fts = Int[]#fts0#Int[]
    #@show faces2delete
    for (i,f) in enumerate(faces0)
      (!in(i,faces2delete)) && (push!(faces,f);push!(fts,fts0[i]))
    end
  else
    faces = faces0
    fts = fts0
  end
  allfaces=vcat(allfaces,faces)
  facetypes = vcat(facetypes,fts)
  depths = vcat(depths,depth*(ones(Int,length(faces))))
  # @show numsides
  # @show numflips
  # @show scale
  # @show centers
  # @show depth
  # @show allfaces
  # @show depths
  # @show facetypes
  for (ctr,scl) in zip(centers,scales)
    a,d,ft = genrecurhex(numsides,numflips,scl,depth-1,ctr; probrecur = probrecur, initial_config = initial_config)#,Curves(),depths,facetypes)#allfaces,depths,facetypes)
    allfaces = vcat(allfaces,a)
    depths = vcat(depths,d)#epth*(ones(Int,length(a))))
    facetypes = vcat(facetypes,ft)
  end
  # @show allfaces
  # @show depths
  # @show facetypes
  #FC,DP,FT =
  #else
  #    cent = (0,0)
  #       genrecurhex(numsides,numsides,(1/numsides)*scl,cent,0,allfaces,depths,facetypes)
  return allfaces,depths,facetypes
end

# function genhex(numsides,numflips,vertvalues,diffs,flippable)
#     for i in 1:numflips
#         fv = rand(map(i->i,flippable))
#         vertvalues,diffs,flippable = flipvert(fv,vertvalues,diffs,flippable,numsides)
#     end
#     faces,facetypes = makefaces(vertvalues,numsides)
#     return faces,facetypes,vertvalues,diffs,flippable
# end
# function genrecursivehex(numsides::Vector{Int},numflips::Int)
#     faces = Vector{Array{Float64}}()
#     vertvalues, diffs, flippable = geninitial(numsides[1])
#     for i in 1:numflips
#         fv = rand(map(i->i,flippable))
#         vertvalues,diffs,flippable = flipvert(fv,vertvalues,diffs,flippable,numsides)
#     end
#     faces,facetypes = makefaces(vertvalues,numsides,scl,center)

# using Gadfly

# function plothex(faces)
#     L = Vector{Gadfly.Layer}()
#     for f in faces
#         l = Gadfly.layer(x=f[:,1],y=f[:,2],Geom.path)[1]
#         push!(L,l)
#     end
#     plot(L)
# end

# numsides = 80
# numflips = 30000
# @time faces,facetypes,vv,dfs,flp,nflp = genhex(numsides, numflips)
