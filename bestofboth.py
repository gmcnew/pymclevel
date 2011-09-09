import math
import optparse
import os
import random
import re
import shutil
import sys
import textwrap
import zlib

import materials
import mclevel
from mclevelbase import ChunkNotPresent

VERSION_STRING = "0.1"

WATER_HEIGHT = 63
MAX_HEIGHT   = 128

class ErosionTask:
    @staticmethod
    def fromString(string):
        (erosionType, direction, posX, posZ) = re.split("\s+", string)
        if erosionType == "corner":
            return CornerErosionTask(Erode.Map.index(direction), int(posX), int(posZ))
        elif erosionType == "edge":
            return EdgeErosionTask(Erode.Map.index(direction), int(posX), int(posZ))
        else:
            raise Exception("unrecognized erosion type '%s'" % (erosionType))
        #    (chunkX, chunkZ, erodeType, xPos, zPos, xMin, xMax, zMin, zMax) = [int(x) for x in re.split("\s+", line)]
        #      if type
        
    def __init__(self, posX, posZ):
        self.posX = posX
        self.posZ = posZ
    
    # Returns the altitude at which the water column starts (assuming it
    # ends at sea level).
    def waterDepth(self, chunk, x, z):
        h = WATER_HEIGHT
        while chunk.Blocks[x, z, h] in [iceID, waterID]:
            h -= 1
        return h
    
    def getChunksAndWaterDepth(self, level):
        chunksToEdit = {}
        chunkWaterDepths = []
        
        for cx in range(-1, 1):
            for cz in range(-1, 1):
                chunk = level.getChunk(self.posX / 16 + cx, self.posZ / 16 + cz)
                chunksToEdit[(cx, cz)] = chunk
                # Find the average water depth along the two borders of this
                # chunk that contact the other chunks involved in this
                # erosion task. For example, if this is the top-left chunk,
                # find the average water depth along its right and bottom
                # edges.
                rowX = 0 if cx == 0 else 15
                rowZ = 0 if cz == 0 else 15
                sumWaterDepths = 0
                for x in range(0, 16):
                    sumWaterDepths += self.waterDepth(chunk, x, rowZ)
                for z in range(0, 16):
                    sumWaterDepths += self.waterDepth(chunk, rowX, z)
            
                chunkWaterDepths.append(sumWaterDepths / 32)
        
        deepestWaterDepth = WATER_HEIGHT
        for wd in chunkWaterDepths:
            if wd < deepestWaterDepth:
                deepestWaterDepth = wd
        
        return (chunksToEdit, deepestWaterDepth)
    
    def run(self, level, erosionWidth, waterWidth):
        raise Exception("not implemented")
    
    # Erodes a column of terrain 1 block wide and 1 block long.
    #
    # relativeDistance: ranges from -1 (on the far edge of the erosion
    #           area) to 1 (on the near edge of the erosion area), with
    #           0 meaning the point is in the middle of the river
    # waterWidth: the width of the river, relative to the width of the erosion
    #           area
    def erode(self, chunk, x, z, relativeDistance, waterWidth, deepestWaterDepth):
        #print("setting (%d,%d) to h=%d" % (x, z, h))
        
        if deepestWaterDepth < WATER_HEIGHT and abs(relativeDistance) < waterWidth:
            # We're in the water and need to slope downward to the ocean
            # floor at deepestWaterDepth
            
            # relativeDistance is on the interval [-waterWidth..waterWidth]
            # (with -waterWidth corresponding to the ocean side and waterWidth
            # corresponding to the ocean side).
            # distanceAcross will be on the interval [0..1] (with 0 being the
            # near side and 1 being the ocean side).
            distanceAcross = (relativeDistance - waterWidth) / (-2 * waterWidth)
            
            # h should range from WATER_HEIGHT to deepestWaterDepth, depending on
            # how far across the river we are.
            h = (deepestWaterDepth - WATER_HEIGHT) * distanceAcross
            h += WATER_HEIGHT
            
        else:
            # We're on land. Make relativeDistance positive for simplicity.
            if relativeDistance < 0:
                relativeDistance *= -1
            relativeRiverDistance = (relativeDistance - waterWidth) / (1 - relativeDistance)
            
            # relativeRiverDistance measures the relative distance from the
            # river's edge to the high point (instead of from the river's center
            # to the high point).
            
            currentTerrainHeight = MAX_HEIGHT - 1
            while chunk.Blocks[x, z, currentTerrainHeight] in leafAndAirIDs:
                currentTerrainHeight -= 1
            
            # relativeDistance is on the interval [0..1],
            # so h will be on the interval [2^0..2^1], which is [1..2].
            #h = 2 ** relativeDistance
            
            # Shift h to the interval [0..1].
            #h -= 1
            #h *= (currentTerrainHeight - WATER_HEIGHT)
            
            h = (currentTerrainHeight - WATER_HEIGHT) * relativeRiverDistance
            if (h < 0):
                h = 0
            h += WATER_HEIGHT
        
        h = int(h)

        #print("(%d,%d) %d" % (x, z, h))
        
        chunkChanged = False
        
        blockWasIce = (chunk.Blocks[x, z, WATER_HEIGHT] == iceID)
        
        if h < MAX_HEIGHT:
        
            # The height at which air will begin.
            airHeight = max(h, WATER_HEIGHT + 1)
            
            # If a tree is standing on terrain that will be preserved,
            # preserve the tree, too.
            while chunk.Blocks[x, z, airHeight] in logIDs:
                airHeight += 1
            
            surfaceHeight = 127
            # We need to dig through the snow layer, since it might just be on
            # top of some leaves.
            while chunk.Blocks[x, z, surfaceHeight] in (leafAndAirIDs + logIDs + [snowLayerID]):
                surfaceHeight -= 1
            
            # Move back up one block if we're under a snow layer.
            if chunk.Blocks[x, z, surfaceHeight + 1] == snowLayerID:
                surfaceHeight += 1
            
            # surfaceHeight is now at ground level. If a tree is standing
            # here, increment surfaceHeight so that the bottommost log will
            # be part of the shifted column (and we can try planting a sapling).
            if chunk.Blocks[x, z, surfaceHeight + 1] in logIDs:
                surfaceHeight += 1
            
            # surfaceHeight is the height of the air block immediately
            # above the surface of the world.
            surfaceHeight += 1
            
            # The terrain between airHeight and surfaceHeight needs to be turned
            # into air.
            
            # Slide terrain downward until it's below airHeight.
            removeDepth = surfaceHeight - airHeight
            removeDepth = min(removeDepth, (airHeight - WATER_HEIGHT) - 1)
            
            # Slide a slice of terrain downward, starting with the bottom-most
            # block.
            for y in range(-removeDepth, 0):
                newBlockID = chunk.Blocks[x, z, surfaceHeight + y]
                
                # The block beneath this block's new location.
                blockBeneathID = chunk.Blocks[x, z, (airHeight + y) - 1]
                
                # If we're moving a snow layer, make sure the block beneath it
                # will be solid. (Some blocks are turned to air during erosion.)
                if newBlockID == snowLayerID:
                    if blockBeneathID == airID:
                        newBlockID = airID
                
                # Don't shift trees downward with the rest of the terrain.
                # This would be difficult to do, since adjacent 1x1 columns are
                # almost always shifted by different amounts. Instead, we'll try
                # to be responsible citizens and plant a new tree for each one
                # that we kill during the erosion process.
                elif newBlockID in logIDs:
                    if blockBeneathID == dirtID:
                        # This log was on dirt. Replace it with a sapling of the
                        # correct type.
                        newBlockID = sapling[newBlockID]
                        print("using a sapling!")
                    else:
                        # This log wasn't on dirt. Turn it into air.
                        newBlockID = airID
                
                # Don't shift leaves, since that could cause weird-looking
                # trees. (If the leaves aren't next to a nearby log after
                # the shift, they'll just decay anyway.)
                elif newBlockID in leafIDs:
                    newBlockID = airID
                
                chunk.Blocks[x, z, airHeight + y] = newBlockID
                chunkChanged = True
            
            # Turn everything in this vertical column into air, but
            # leave leaves alone (to avoid weird-looking half-trees). Logs may
            # be destroyed during erosion, so some leaves may no longer be near
            # any logs, but these leaves will decay naturally in the game.
            for ah in range(airHeight, 128):
                if chunk.Blocks[x, z, ah] not in leafIDs:
                    chunk.Blocks[x, z, ah] = airID
                    chunkChanged = True
            if h <= WATER_HEIGHT + 1:
                if h <= WATER_HEIGHT:
                    chunk.Blocks[x, z, h : WATER_HEIGHT + 1] = waterID
                    chunkChanged = True
                # Turn non-water, non-ice blocks along the shoreline, or under the water, into sand.
                if chunk.Blocks[x, z, h - 1] != waterID \
                        and chunk.Blocks[x, z, h - 1] != iceID:
                    chunk.Blocks[x, z, h - 1] = sandID
                    chunkChanged = True
        
        if blockWasIce:
            # Restore ice that was far from the center of the river.
            # A larger relative distance from the center of the river should
            # result in a greater chance of restoring the block of ice.
            if random.random() < abs(relativeDistance):
                chunk.Blocks[x, z, WATER_HEIGHT] = iceID
        
        return chunkChanged
        
class CornerErosionTask(ErosionTask):
    def __init__(self, cornerDirection, cornerPosX, cornerPosZ):
        ErosionTask.__init__(self, cornerPosX, cornerPosZ)
        self.cornerDirection = cornerDirection
        
    def __repr__(self):
        return "corner %-2s %d %d" % (Erode.Map[self.cornerDirection], self.posX, self.posZ)
        
    def run(self, level, erosionWidth = 8, waterWidth = 3):
        try:
            (chunksToEdit, deepestWaterDepth) = self.getChunksAndWaterDepth(level)
        except ChunkNotPresent:
            return False
            
        if self.cornerDirection == Erode.TL:
            highPoint = (8, 8)
        elif self.cornerDirection == Erode.TR:
            highPoint = (-8, 8)
        elif self.cornerDirection == Erode.BR:
            highPoint = (-8, -8)
        elif self.cornerDirection == Erode.BL:
            highPoint = (8, -8)
            
        for cx in range(-1, 1):
            for cz in range(-1, 1):
                chunkChanged = False
                chunk = chunksToEdit[(cx, cz)]
                highPointX = highPoint[0] - (cx * 16)
                highPointZ = highPoint[1] - (cz * 16)
                for x in range(-8 * cx, -8 * (cx - 1)):
                    for z in range(-8 * cz, -8 * (cz - 1)):
                        dx = x - (highPointX - 0.5)
                        dz = z - (highPointZ - 0.5)
                        
                        distanceFromCenter = math.sqrt(dx * dx + dz * dz)
                        distanceFromEdge = 8 - distanceFromCenter
                        
                        if abs(distanceFromEdge) < erosionWidth:
                            relativeDistance = float(distanceFromEdge) / float(erosionWidth)
                            chunkChanged |= self.erode(chunk, x, z, relativeDistance, .375, deepestWaterDepth)

                if chunkChanged:
                    chunk.chunkChanged()
        return True
        
class EdgeErosionTask(ErosionTask):
    def __init__(self, edgeDirection, edgePosX, edgePosZ):
        ErosionTask.__init__(self, edgePosX, edgePosZ)
        self.edgeDirection = edgeDirection
        
    def __repr__(self):
        return "edge   %-2s %d %d" % (Erode.Map[self.edgeDirection], self.posX, self.posZ)
        
    def run(self, level, erosionWidth = 8, waterWidth = 3):
        try:
            (chunksToEdit, deepestWaterDepth) = self.getChunksAndWaterDepth(level)
        except ChunkNotPresent:
            return False
        
        for cx in range(-1, 1):
            for cz in range(-1, 1):
                chunkChanged = False
                xMin = -8 * cx
                xMax = -8 * (cx - 1)
                zMin = -8 * cz
                zMax = -8 * (cz - 1)
                chunk = chunksToEdit[(cx, cz)]
                for x in range(xMin, xMax):
                    for z in range(zMin, zMax):
                        if self.edgeDirection == Erode.HE:
                            # horizontal edge
                            distanceFromCenter = abs(7.5 - z)
                        elif self.edgeDirection == Erode.VE:
                            # vertical edge
                            distanceFromCenter = abs(7.5 - x)
                        else:
                            raise Exception("unrecognized edge direction %d (%s)" % (self.edgeDirection, Erode.Map[self.edgeDirection]))
                            
                        distanceFromEdge = 8 - distanceFromCenter
                        
                        if abs(distanceFromEdge) < erosionWidth:
                            relativeDistance = float(distanceFromEdge) / float(erosionWidth)
                            chunkChanged |= self.erode(chunk, x, z, relativeDistance, .375, deepestWaterDepth)
                            
                if chunkChanged:
                    chunk.chunkChanged()
        return True


# These values indicate the shape a chunk should have after erosion.
# T, B, L, and R mean "top", "bottom", "left" and "right", respectively.
# TL, TR, BL, and BR refer to the four kinds of corners: "top left",
# "top right", "bottom left", and "bottom right", respectively.
class Erode:
    Map = [ "TL", "TR", "BL", "BR", "VE", "HE" ]
    TL  = 0 # top-left corner
    TR  = 1 # top-right corner
    BL  = 2 # bottom-left corner
    BR  = 3 # bottom-right corner
    VE  = 4 # vertical edge
    HE  = 5 # horizontal edge

airID       = materials.materials.Air.ID
dirtID      = materials.materials.Dirt.ID
iceID       = materials.materials.Ice.ID
sandID      = materials.materials.Sand.ID
snowLayerID = materials.materials.SnowLayer.ID
waterID     = materials.materials.WaterStill.ID
    
logIDs = [  
    materials.materials.Wood.ID,
    materials.materials.Ironwood.ID,
    materials.materials.BirchWood.ID,
    ]

# Map log IDs to the corresponding sapling types. This is used when replanting
# trees.
sapling = {
    materials.materials.Wood.ID:        materials.materials.Sapling.ID,
    materials.materials.Ironwood.ID:    materials.materials.SpruceSapling.ID,
    materials.materials.BirchWood.ID:   materials.materials.BirchSapling.ID,
}

leafIDs = [
    materials.materials.Leaves.ID,
    materials.materials.PineLeaves.ID,
    materials.materials.BirchLeaves.ID,
    ]

leafAndAirIDs = [
    materials.materials.Air.ID,
    materials.materials.Leaves.ID,
    materials.materials.PineLeaves.ID,
    materials.materials.BirchLeaves.ID,
    ]

def find_edges(worldDir, edgeFilename):
    level = mclevel.fromFile(worldDir)
    edgeFile = open(edgeFilename, "w")
    sys.stdout.write("finding edges...")
    
    chunks = []
    
    for chunk in level.allChunks:
        chunks.append(chunk)
    
    erodeTasks = []
    
    examined = 0
    lastProgress = 0
    numChunks = len(chunks)
    
    for chunk in chunks:
        checkChunk(level, chunk, erodeTasks)
        examined += 1
        progress = examined * 100 / numChunks
        if progress != lastProgress:
            lastProgress = progress
            sys.stdout.write("\rfinding edges (%d%%)..." % (progress))
    print("")
    
    edgeFile.write("# erodeType erodeDirection posX posZ\n")
    
    numEdgeChunks = 0
    for task in erodeTasks:
        edgeFile.write("%s\n" % (task))
        numEdgeChunks += 1
    edgeFile.close()
    print("found %d edge(s)" % (numEdgeChunks))

def smooth(worldDir, edgeFilename, width = 16):
    level = mclevel.fromFile(worldDir)
    newEdgeFile = open(edgeFilename + ".tmp", "w")
    edgeFile = open(edgeFilename, "r")
    
    width = int(width) / 2
    
    erosionTasks = []
    
    for line in edgeFile.readlines():
        originalLine = line
        line = line.strip()
        # Preserve comments
        if line.startswith("#"):
            newEdgeFile.write(originalLine)
        else:
            task = ErosionTask.fromString(line)
            erosionTasks.append(task)
    
    edgeFile.close()
    
    numTasks = len(erosionTasks)
    
    skipped = 0
    smoothed = 0
    
    if erosionTasks:
        examined = 0
        for erosionTask in erosionTasks:
            examined += 1
            sys.stdout.write("\rexamining edge %d of %d..." % (examined, numTasks))
            
            # If the task didn't run (because it requires chunks that
            # haven't been generated yet), write it back to edges.txt.
            if erosionTask.run(level, width):
                smoothed += 1
            else:
                skipped += 1
                newEdgeFile.write("%s\n" % (task))
            
        
        level.saveInPlace()
        
        print("")
    
    newEdgeFile.close()
    
    if smoothed:
        print("smoothed %d edge(s)" % (smoothed))
        shutil.move(newEdgeFile.name, edgeFilename)
    else:
        os.remove(newEdgeFile.name)
    
    if skipped:
        print("%d edge(s) can't be smoothed yet, since they're not fully explored" % (skipped))
    elif smoothed == numTasks:
        print("the map is perfectly smoothed -- nothing to do!")

def addCorner(chunkPos, erodeList, erodeType):
    (chunkX, chunkZ) = chunkPos
    erodeList.append(
        CornerErosionTask(
            erodeType,
            chunkX * 16,
            chunkZ * 16
        )
    )

def addEdge(chunkPos, erodeList, erodeType):
    (chunkX, chunkZ) = chunkPos
    erodeList.append(
        EdgeErosionTask(
            erodeType,
            chunkX * 16,
            chunkZ * 16
        )
    )

# Examine a chunk in a level. For each edge that's found, add a
# (chunk, direction) pair to erodeQueue. Return the number of pairs
# added to the queue.
def checkChunk(level, coords, toErode):
    
    aroundMe = [(-1, -1), (0, -1), (1, -1),
                (-1,  0),          (1,  0),
                (-1,  1), (0,  1), (1,  1)]
    
    (TL, T, TR, L, R, BL, B, BR) = range(0, 8)
    
    neighbors = [True] * 8
    
    onPerimeter = False
    
    for i in range(len(aroundMe)):
        a = aroundMe[i]
        if (coords[0] + a[0], coords[1] + a[1]) not in level.allChunks:
            onPerimeter = True
            neighbors[i] = False
    
    if onPerimeter:
    
        # Top-left corner
        if not (neighbors[TL] or neighbors[T] or neighbors[L]):
            addCorner(coords, toErode, Erode.TL)
        
        # Top-right corner
        if not (neighbors[T] or neighbors[TR] or neighbors[R]):
            coordsRight = (coords[0] + 1, coords[1])
            addCorner(coordsRight, toErode, Erode.TR)
        
        # Bottom-right corner
        if not (neighbors[R] or neighbors[BR] or neighbors[B]):
            coordsBelowAndRight = (coords[0] + 1, coords[1] + 1)
            addCorner(coordsBelowAndRight, toErode, Erode.BR)
        
        # Bottom-left corner
        if not (neighbors[B] or neighbors[BL] or neighbors[L]):
            coordsBelow = (coords[0], coords[1] + 1)
            addCorner(coordsBelow, toErode, Erode.BL)
        
        # Top-left corner (inverted)
        if neighbors[L] and neighbors[T] and not neighbors[TL]:
            addCorner(coords, toErode, Erode.BR)
        
        # Top-right corner (inverted)
        if neighbors[T] and neighbors[R] and not neighbors[TR]:
            coordsRight = (coords[0] + 1, coords[1])
            addCorner(coordsRight, toErode, Erode.BL)
        
        # Bottom-right corner (inverted)
        if neighbors[R] and neighbors[B] and not neighbors[BR]:
            coordsBelowAndRight = (coords[0] + 1, coords[1] + 1)
            addCorner(coordsBelowAndRight, toErode, Erode.TL)
        
        # Bottom-left corner (inverted)
        if neighbors[B] and neighbors[L] and not neighbors[BL]:
            coordsBelow = (coords[0], coords[1] + 1)
            addCorner(coordsBelow, toErode, Erode.TR)
        
        if neighbors[L]:
            if not (neighbors[T] or neighbors[TL]):
                addEdge(coords, toErode, Erode.HE)
            
            if not (neighbors[B] or neighbors[BL]):
                coordsBelow = (coords[0], coords[1] + 1)
                addEdge(coordsBelow, toErode, Erode.HE)
        
        if neighbors[T]:
            if not (neighbors[L] or neighbors[TL]):
                addEdge(coords, toErode, Erode.VE)
        
            if not (neighbors[R] or neighbors[TR]):
                coordsRight = (coords[0] + 1, coords[1])
                addEdge(coordsRight, toErode, Erode.VE)
                
        # These "checkerboard" cases are kind of weird.
        if neighbors[TL] and not (neighbors[L] or neighbors[T]):
            addEdge(coords, toErode, Erode.HE)
            addEdge(coords, toErode, Erode.VE)
        
        if neighbors[TR] and not (neighbors[R] or neighbors[T]):
            coordsRight = (coords[0] + 1, coords[1])
            addEdge(coordsRight, toErode, Erode.HE)
            addEdge(coordsRight, toErode, Erode.VE)
    
    return len(toErode)

def get_info_text():
    return "\n".join([
        "bestofboth version %s" % VERSION_STRING,
        "(a tool for smoothing terrain discontinuities in Minecraft worlds)",
        "http://github.com/gmcnew/pymclevel",
        ""])

def get_usage_text():
    usage = """
    bestofboth --find-edges <path_to_world>
    bestofboth --smooth <path_to_world>"""

    usageWithSmooth = """
    bestofboth --find-edges <path_to_world>
    bestofboth --smooth <path_to_world> [--width <1-16>]"""
    
    # A paragraph is a list of lines.
    paragraphs = [[usage]]
    
    paragraphs.append(textwrap.wrap(
        "This script must be run in two steps. The first is the " \
        "--find-edges step, which examines a world and finds its edges. " \
        "Next is the --smooth step, which smooths edges by carving a river " \
        "between old chunks and newly-generated ones."))
    
    paragraphs.append(textwrap.wrap(
        "You can run the --smooth step multiple times as players explore " \
        "edges and cause new chunks to be generated along them. Eventually, " \
        "if all chunks next to edges have been generated, the script will " \
        "report that the map is perfectly smoothed. At this point further use " \
        "of the script is unnecessary."))
    
    paragraphs.append(
        ["Typical use:"] +
        [("    %s" % x) for x in [
            "bestofboth --find-edges <path_to_world>",
            "<upgrade Minecraft to a version with new terrain generation code>",
            "...",
            "<play Minecraft, explore, and cause new terrain to be generated>",
            "bestofboth --smooth <path_to_world>",
            "...",
            "<more exploration and edge discovery>",
            "bestofboth --smooth <path_to_world>",
            "...",
            "<finish exploring edges and new terrain along them>",
            "bestofboth --smooth <path_to_world>",
            ]
        ]
    )
    
    return "\n\n".join(["\n".join(p) for p in paragraphs])
    
def main():
    
    random.seed(0)
    
    parser = optparse.OptionParser(usage = get_usage_text())
    parser.add_option("--find-edges",
                    dest="find_edges",
                    metavar = "path",
                    help="path to the world to examine")
    parser.add_option("--smooth",
                    dest="smooth",
                    metavar = "path",
                    help="path to the world to smooth")
    """
    parser.add_option("--width", dest="width", 
                    default = "16",
                    help="width of the river")
    """
    
    print(get_info_text())

    (options, args) = parser.parse_args()
    
    worldDir = options.find_edges or options.smooth
    if worldDir:
        edgeFilePath = os.path.join(worldDir, "edges.txt")
    
    errorText = None
    
    if options.find_edges and options.smooth:
        errorText = "--find-edges and --smooth can't be specified " \
            "at the same time. Please run with --find-edges first, " \
            "then run with --smooth."
    elif not (options.find_edges or options.smooth):
        parser.print_help()
    elif not os.path.exists(os.path.join(worldDir, "level.dat")):
        errorText = "'%s' is not a Minecraft world directory (no " \
            "level.dat file was found)." % (worldDir)
    elif options.smooth and not os.path.exists(edgeFilePath):
        errorText = "Edge file '%s' does not exist. Run with " \
            "--find-edges to create the edge file, which must exist when " \
            "--smooth is specified." \
            % (edgeFilePath)
    elif options.find_edges and os.path.exists(edgeFilePath):
        errorText = "Edge file '%s' already exists. Did you mean " \
            "to specify --smooth instead?" \
            % (edgeFilePath)
    
    
    if errorText:
        parser.error("\n" + "\n".join(textwrap.wrap(errorText)))
    
    """
    elif options.width and (int(options.width) < 1 or int(options.width) > 16):
        parser.error("--width must be between 1 and 16 (inclusive)")
    """
    
    # Phew! Now that the arguments have been validated...
    if options.find_edges:
        find_edges(worldDir, edgeFilePath)
    elif options.smooth:
        # TODO: Fix the "--width" argument.
        #smooth(options.smooth, edgeFilePath, options.width)
        smooth(worldDir, edgeFilePath)

if __name__ == "__main__":
    main()
