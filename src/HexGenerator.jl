module HexGenerator

using PointsNCurves

export genhex,
       geninitial0,
       geninitial,
       geninitial_alternate,
       flipvert,
       makefaces,
       makeface,
       ind2coords,
       vertneighbors,
       faceneighbors,
       face2verts,
       facecoords,
       genrecurhex,
       numvertsthisrow,
       numfacesthisrow,
       center,
       compute_diffs_and_flippable

include("makehexes.jl")
end
