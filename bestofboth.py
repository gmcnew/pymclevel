import math
import optparse
import os
import Queue
import random
import re
import shutil
import tempfile

import materials
import mclevel
from mclevelbase import ChunkNotPresent

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
            relativeDistance = (relativeDistance - waterWidth) / (1 - relativeDistance)
            
            # relativeDistance now measures the relative distance from the
            # river's edge to the high point (instead of from the river's center
            # to the high point).
            
            currentTerrainHeight = MAX_HEIGHT - 1
            while chunk.Blocks[x, z, currentTerrainHeight] in leafIDs \
                    or chunk.Blocks[x, z, currentTerrainHeight] == airID:
                currentTerrainHeight -= 1
            
            # relativeDistance is on the interval [0..1],
            # so h will be on the interval [2^0..2^1], which is [1..2].
            #h = 2 ** relativeDistance
            
            # Shift h to the interval [0..1].
            #h -= 1
            #h *= (currentTerrainHeight - WATER_HEIGHT)
            
            h = (currentTerrainHeight - WATER_HEIGHT) * relativeDistance
            if (h < 0):
                h = 0
            h += WATER_HEIGHT
        
        h = int(h)

        #print("(%d,%d) %d" % (x, z, h))
        
        chunkChanged = False
        
        if h < MAX_HEIGHT:
        
            # The height at which air will begin.
            ah = h
            
            # If a tree is standing on terrain that will be preserved,
            # preserve the tree, too.
            while chunk.Blocks[x, z, ah] in logIDs:
                ah += 1
            
            # Turn everything in this vertical column into air, but
            # leave leaves alone (to avoid weird-looking half-trees).
            for h2 in range(ah, 128):
                if chunk.Blocks[x, z, h2] not in leafIDs:
                    chunk.Blocks[x, z, h2] = airID
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
# T, B, L and R mean "top", "bottom", "left" and "right", respectively.
# There are two kinds of corners, illustrated by this fine piece of
# artwork:
#
#   a1 a2
#   a3 b4
#
# Imagine maps 'a' and 'b' are being merged. The top-left corner of 'b4'
# needs to be eroded. This corresponds to "TL" below. The bottom-right
# corner of 'a1' needs to be eroded as well. However, the shape of this
# erosion is inverted -- 'a1' is a concave point on 'a', not a convex
# point. This corresponds to "NBR" below -- the "N" meaning "inverted"
# and "BL" meaning "bottom-right".
class Erode:
    Map = [ "TL", "TR", "BL", "BR", "VE", "HE" ]
    TL  = 0 # top-left corner
    TR  = 1 # top-right corner
    BL  = 2 # bottom-left corner
    BR  = 3 # bottom-right corner
    VE  = 4 # vertical edge
    HE  = 5 # horizontal edge
TL = 0
T  = 1
TR = 2
L  = 3
R  = 4
BL = 5
B  = 6
BR = 7

airID   = materials.materials.Air.ID
iceID   = materials.materials.Ice.ID
sandID  = materials.materials.Sand.ID
waterID = materials.materials.WaterStill.ID
    
logIDs = [  
    materials.materials.Wood.ID,
    materials.materials.Ironwood.ID,
    materials.materials.BirchWood.ID,
    ]

leafIDs = [
    materials.materials.Leaves.ID,
    materials.materials.PineLeaves.ID,
    materials.materials.BirchLeaves.ID,
    ]

def find_edges(worldDir, edgeFilename):
    level = mclevel.fromFile(worldDir)
    edgeFile = open(edgeFilename, "w")
    print("world: %s" % (worldDir))
    print("finding edges...")
    
    erodeQ = Queue.Queue()
    
    for chunk in level.allChunks:
        checkChunk(level, chunk, erodeQ)
    
    edgeFile.write("# erodeType erodeDirection posX posZ")
    
    numEdgeChunks = 0
    while not erodeQ.empty():
        task = erodeQ.get(True)
        edgeFile.write("%s\n" % (task))
        numEdgeChunks += 1
    edgeFile.close()
    print("found %d perimeter features" % (numEdgeChunks))

def smooth(worldDir, width, edgeFilename):
    level = mclevel.fromFile(worldDir)
    newEdgeFile = open(edgeFilename + ".tmp", "w")
    edgeFile = open(edgeFilename, "r")
    
    width = int(width) / 2
    
    erosionTasks = []
    
    numChunks = 0
    for line in edgeFile.readlines():
        originalLine = line
        line = line.strip()
        # Preserve comments
        if line.startswith("#"):
            newEdgeFile.write(originalLine)
        else:
            task = ErosionTask.fromString(line)
            numChunks += 1
            erosionTasks.append(task)
    
    edgeFile.close()
    
    numEroded = len(erosionTasks)
    
    if erosionTasks:
        print("smoothing %d perimeter feature(s)..." % (numEroded))
        
        i = 0
        lastReportedProgress = 0
        for erosionTask in erosionTasks:
            i += 1
            progress = i * 100 / numEroded
            if progress >= lastReportedProgress + 10:
                lastReportedProgress = progress
                print("%d%%" % (progress))
            
            # If the task didn't run (because it requires chunks that
            # haven't been generated yet), write it back to edges.txt.
            if not erosionTask.run(level, width):
                newEdgeFile.write("%s\n" % (task))
            
        
        level.saveInPlace()
        print("smoothed %d perimeter feature(s)" % (numEroded))
    
    newEdgeFile.close()
    
    if numChunks - numEroded:
        print("%d perimeter feature(s) were not smoothed (since they haven't been explored completely)" % (numChunks - numEroded))
    elif not numChunks:
        print("the map is perfectly smoothed -- nothing to do!")
    
    shutil.move(newEdgeFile.name, edgeFilename)

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
def checkChunk(level, coords, erodeQueue):
    
    tasksAdded = 0
    
    aroundMe = [(-1, -1), (0, -1), (1, -1),
                (-1,  0),          (1,  0),
                (-1,  1), (0,  1), (1,  1)]
    
    neighbors = [True] * 8
    
    onPerimeter = False
    
    for i in range(len(aroundMe)):
        a = aroundMe[i]
        if (coords[0] + a[0], coords[1] + a[1]) not in level.allChunks:
            onPerimeter = True
            neighbors[i] = False
    
    toErode = []
    
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
        
    for task in toErode:
        erodeQueue.put(task)
        tasksAdded += 1
    
    return tasksAdded

def main():
    usage = """
bestofboth --find-edges <path_to_world>
bestofboth --smooth <path_to_world> [--width <1-16>]
"""
    parser = optparse.OptionParser(usage = usage)
    parser.add_option("--find-edges", dest="find_edges",
                    help="world to examine")
    parser.add_option("--smooth", dest="smooth", 
                    help="world to smooth")
    parser.add_option("--width", dest="width", 
                    default = "16",
                    help="width of the river")

    (options, args) = parser.parse_args()
    
    worldDir = options.find_edges or options.smooth
    
    if options.find_edges and options.smooth:
        parser.error("--find-edges and --smooth can't be specified " \
            "at the same time. Please run with --find-edges first, " \
            "then run with --smooth.")
    elif not (options.find_edges or options.smooth):
        parser.error("Must specify --find-edges or --smooth.")
    elif not os.path.exists(os.path.join(worldDir, "level.dat")):
            parser.error("'%s' is not a Minecraft world directory (no " \
                "level.dat file was found)." % (options.worldDir))
    elif options.smooth and not os.path.exists(os.path.join(options.smooth, "edges.txt")):
        parser.error("Edge file '%s' does not exist. Run with " \
            "--find-edges first. The edge file must exist when " \
            "--smooth is specified." \
            % (os.path.join(options.smooth, "edges.txt")))
    elif options.width and (int(options.width) < 1 or int(options.width) > 16):
        parser.error("--width must be between 1 and 16 (inclusive)")
    
    edgeFilePath = os.path.join(worldDir, "edges.txt")
    
    # Phew! Now that the arguments have been validated...
    if options.find_edges:
        find_edges(options.find_edges, edgeFilePath)
    elif options.smooth:
        smooth(options.smooth, options.width, edgeFilePath)

if __name__ == "__main__":
    main()
