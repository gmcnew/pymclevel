import math
from optparse import OptionParser
import os
import Queue
import random
import re
import shutil
import sys
import tempfile
import threading
import time

import materials
import mclevel
from mclevelbase import ChunkNotPresent

X = 0
Y = 1
Z = 2

WATER_HEIGHT = 63
MAX_HEIGHT   = 128

NUM_WORKERS = 2


import signal
import sys


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
    Map = [ "TL", "T", "TR", "L", "R", "BL", "B", "BR",
            "NTL", "NTR", "NBL", "NBR" ]
    TL  = 0
    T   = 1
    TR  = 2
    L   = 3 # toward the left
    R   = 4 # toward the right
    BL  = 5 # toward the bottom-left 
    B   = 6
    BR  = 7
    NTL = 8
    NTR = 9
    NBL = 10
    NBR = 11

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


# An Eroder takes (chunk, direction) pairs from a queue (erodeQueue),
# erodes them, and returns them to a different queue (progressQueue).
# However, map smoothing is IO-bound, not CPU-bound, so multithreading
# this probably isn't worth the complication.
class Eroder(threading.Thread):
    
    def __init__(self, ID, level, erodeQueue, progressQueue):
        threading.Thread.__init__(self)
        self.ID = ID
        self.level = level
        self.erodeQueue = erodeQueue
        self.progressQueue = progressQueue
        
    def run(self):
        while True:
            # Wait forever for a new task. Null tasks (where chunk is
            # None) will be at the end of the queue to signal the end of
            # work.
            (coords, direction) = self.erodeQueue.get(True)
            if coords is None:
                print("worker %d is returning" % self.ID)
                break
            if random.random() > 0.999:
                try:
                    chunk = self.level.getChunk(coords[0], coords[1])
                    erode(chunk, direction)
                except ChunkNotPresent as cnp:
                    pass
            self.progressQueue.put(coords)

def find_edges(worldDir, edgeFilename):
    level = mclevel.fromFile(worldDir)
    edgeFile = open(edgeFilename, "w")
    print("world: %s" % (worldDir))
    print("finding edges...")
    
    erodeQ = Queue.Queue()
    
    for chunk in level.allChunks:
        checkChunk(level, chunk, erodeQ)
    
    edgeFile.write("# chunkCoordX chunkCoordY erosionType\n")
    
    numEdgeChunks = 0
    while not erodeQ.empty():
        (coords, direction) = erodeQ.get(True)
        edgeFile.write("%d %d %s\n" % (coords[0], coords[1], Erode.Map[direction]))
        numEdgeChunks += 1
    edgeFile.close()
    print("found %d edge chunk(s)" % (numEdgeChunks))

def smooth(worldDir, width, edgeFilename):
    level = mclevel.fromFile(worldDir)
    edgeFile = open(edgeFilename, "r")
    newEdgeFile = open(edgeFilename + ".tmp", "w")
    
    width = int(width) / 2
    
    toErode = []
    
    numChunks = 0
    for line in edgeFile.readlines():
        originalLine = line
        line = line.strip()
        # Preserve comments
        if line.startswith("#"):
            newEdgeFile.write(originalLine)
        else:
            (chunkX, chunkZ, directionString) = re.split("\s+", line)
            (chunkX, chunkZ) = (int(chunkX), int(chunkZ))
            direction = Erode.Map.index(directionString)
            numChunks += 1
            if direction < 0:
                raise Exception("unrecognized direction '%s' in edge file" % (directionString))
            try:
                chunk = level.getChunk(chunkX, chunkZ)
                toErode.append((chunk, direction))
            except ChunkNotPresent as cnp:
                newEdgeFile.write("%d %d %s\n" % (chunkX, chunkZ, directionString))
    
    newEdgeFile.close()
    edgeFile.close()
    
    numEroded = len(toErode)
    
    if toErode:
        print("smoothing %d chunk(s)..." % (numEroded))
        
        i = 0
        lastReportedProgress = 0
        for (chunk, direction) in toErode:
            i += 1
            progress = i * 100 / numEroded
            if progress >= lastReportedProgress + 10:
                lastReportedProgress = progress
                print("%d%%" % (progress))
            erode(chunk, direction, width)
            chunk.chunkChanged()
        
        level.saveInPlace()
        print("smoothed %d chunk(s)" % (numEroded))
    
    if numChunks - numEroded:
        print("%d chunk(s) were not smoothed (since they haven't been explored)" % (numChunks - numEroded))
    elif not numChunks:
        print("the map is perfectly smoothed -- nothing to do!")
    
    shutil.move(newEdgeFile.name, edgeFilename)

# Examine a chunk in a level. For each edge that's found, add a
# (chunk, direction) pair to erodeQueue. Return the number of pairs
# added to the queue.
def checkChunk(level, chunk, erodeQueue):
    
    tasksAdded = 0
    
    aroundMe = [(-1, -1), (0, -1), (1, -1),
                (-1,  0),          (1,  0),
                (-1,  1), (0,  1), (1,  1)]
    
    neighbors = [True] * 8
    
    onPerimeter = False
    
    for i in range(len(aroundMe)):
        a = aroundMe[i]
        if (chunk[0] + a[0], chunk[1] + a[1]) not in level.allChunks:
            onPerimeter = True
            neighbors[i] = False
    
    toErode = []
    
    if onPerimeter:
    
        # Top-left corner
        if not (neighbors[TL] or neighbors[T] or neighbors[L]):
            toErode.append((chunk, Erode.TL))
            toErode.append(((chunk[0] - 1, chunk[1] - 1), Erode.NBR))
        
        # Top-right corner
        if not (neighbors[T] or neighbors[TR] or neighbors[R]):
            toErode.append((chunk, Erode.TR))
            toErode.append(((chunk[0] + 1, chunk[1] - 1), Erode.NBL))
        
        # Bottom-right corner
        if not (neighbors[R] or neighbors[BR] or neighbors[B]):
            toErode.append((chunk, Erode.BR))
            toErode.append(((chunk[0] + 1, chunk[1] + 1), Erode.NTL))
        
        # Bottom-left corner
        if not (neighbors[B] or neighbors[BL] or neighbors[L]):
            toErode.append((chunk, Erode.BL))
            toErode.append(((chunk[0] - 1, chunk[1] + 1), Erode.NTR))
            
        if not neighbors[T]:
            toErode.append((chunk, Erode.T))
            toErode.append(((chunk[0], chunk[1] - 1), Erode.B))
            
        if not neighbors[B]:
            toErode.append((chunk, Erode.B))
            toErode.append(((chunk[0], chunk[1] + 1), Erode.T))
        
        if not neighbors[L]:
            toErode.append((chunk, Erode.L))
            toErode.append(((chunk[0] - 1, chunk[1]), Erode.R))
        
        if not neighbors[R]:
            toErode.append((chunk, Erode.R))
            toErode.append(((chunk[0] + 1, chunk[1]), Erode.L))
        
        # Top-left corner (inverted)
        if neighbors[T] and neighbors[L] and not neighbors[TL]:
            toErode.append((chunk, Erode.NTL))
            toErode.append(((chunk[0] - 1, chunk[1] - 1), Erode.BR))
        
        # Top-right corner (inverted)
        if neighbors[T] and neighbors[R] and not neighbors[TR]:
            toErode.append((chunk, Erode.NTR))
            toErode.append(((chunk[0] + 1, chunk[1] - 1), Erode.BL))
        
        # Bottom-right corner (inverted)
        if neighbors[B] and neighbors[R] and not neighbors[BR]:
            toErode.append((chunk, Erode.NBR))
            toErode.append(((chunk[0] + 1, chunk[1] + 1), Erode.TL))
        
        # Bottom-left corner (inverted)
        if neighbors[B] and neighbors[L] and not neighbors[BL]:
            toErode.append((chunk, Erode.NBL))
            toErode.append(((chunk[0] - 1, chunk[1] + 1), Erode.TR))
    
    for (chunkCoords, direction) in toErode:
        erodeQueue.put((chunkCoords, direction))
        tasksAdded += 1
    
    return tasksAdded

# erosionWidth: half of the total width of the erosion area, from the
#   middle of the river to one edge
# riverWidth: half of the width, in blocks, of the river running through
#   the erosion area
# maxPreservedHeight: Within the erosion area, everything above this
#   height will be eroded away, even on the edges of the erosion area.
#   There will be a slope between the edges and the river.
def erosionHeight(x, z, direction, erosionWidth = 8, riverWidth = 1, maxPreservedHeight = 90):
    
    # Corners are eroded based on the location of a high point -or- a
    # low point.
    lowPoint = None
    highPoint = None
    if direction == Erode.TL:
        highPoint = (erosionWidth, erosionWidth)
    elif direction == Erode.TR:
        highPoint = (16 - erosionWidth, erosionWidth)
    elif direction == Erode.BL:
        highPoint = (erosionWidth, 16 - erosionWidth)
    elif direction == Erode.BR:
        highPoint = (16 - erosionWidth, 16 - erosionWidth)
    elif direction == Erode.NTL:
        lowPoint = (0, 0)
    elif direction == Erode.NTR:
        lowPoint = (16, 0)
    elif direction == Erode.NBL:
        lowPoint = (0, 16)
    elif direction == Erode.NBR:
        lowPoint = (16, 16)
    
    if highPoint:
        dx = highPoint[0] - x - 0.5
        dz = highPoint[1] - z - 0.5
        distanceFromCenter = math.sqrt((dx * dx) + (dz * dz))
        distanceFromEdge = erosionWidth - distanceFromCenter
    elif lowPoint:
        dx = lowPoint[0] - x - 0.5
        dz = lowPoint[1] - z - 0.5
        distanceFromEdge = math.sqrt((dx * dx) + (dz * dz))
    else:
        if direction == Erode.T:
            distanceFromEdge = z + 0.5
        elif direction == Erode.B:
            distanceFromEdge = 15.5 - z
        elif direction == Erode.L:
            distanceFromEdge = x + 0.5
        elif direction == Erode.R:
            distanceFromEdge = 15.5 - x
        
    if distanceFromEdge < 0:
        distanceFromEdge = 0
    
    # We're far enough from the edge that no erosion should occur.
    if distanceFromEdge > erosionWidth:
        return MAX_HEIGHT
    
    distanceFromRiver = distanceFromEdge - riverWidth
    if distanceFromRiver < 0:
        distanceFromRiver = 0
    
    # (distanceFromEdge / halfErosionWidth) is on the interval [0..1],
    # so h will be on the interval [2^0..2^1], or [1..2].
    h = 2 ** (float(distanceFromRiver) / float(erosionWidth - riverWidth))
    
    #if direction == Erode.T and z == 1:
        #print(h)
    
    # Shift h to the interval [0..1].
    h -= 1
    #if direction == Erode.T and z == 1:
        #print(h)
    
    h *= maxPreservedHeight - WATER_HEIGHT
    #if direction == Erode.T and z == 1:
        #print(h)
    
    h += WATER_HEIGHT
    #if direction == Erode.T and z == 1:
        #print(h)
    
    #if direction == Erode.T and z == 1:
        #print("%d,%d: %d %f %f %f" % (x, z, h, distanceFromRiver, distance, halfRiverWidth))
        #print("")
    
    """
    if (x == 1 and z == 7) or (x == 14 and z == 7) or (x == 1 and z == 8) or (x == 14 and z == 8):
        if direction in [Erode.TL, Erode.NTL, Erode.BR, Erode.NBR, Erode.TR, Erode.NTR, Erode.BL, Erode.NBL]:
            print("%s %d,%d: %f" % (Erode.Map[direction], x, z, h))
    """
    
    """
    if x == 0 and z == 0:
        if direction == Erode.TR:
            print("TR %d,%d: %f" % (x, z, h))
        elif direction == Erode.TL:
            print("TL %d,%d: %f" % (x, z, h))
        elif direction == Erode.BL:
            print("BL %d,%d: %f" % (x, z, h))
        elif direction == Erode.BR:
            print("BR %d,%d: %f" % (x, z, h))
    """
    
    return int(h)

def erosionBoundaries(direction, distance):
    if direction == Erode.TL or direction == Erode.NTL:
        return (0, distance, 0, distance)
    elif direction == Erode.TR or direction == Erode.NTR:
        return (16 - distance, 16, 0, distance)
    elif direction == Erode.BR or direction == Erode.NBR:
        return (16 - distance, 16, 16 - distance, 16)
    elif direction == Erode.BL or direction == Erode.NBL:
        return (0, distance, 16 - distance, 16)
    elif direction == Erode.T:
        return (0, 16, 0, distance)
    elif direction == Erode.B:
        return (0, 16, 16 - distance, 16)
    elif direction == Erode.L:
        return (0, distance, 0, 16)
    elif direction == Erode.R:
        return (16 - distance, 16, 0, 16)

def erode(chunk, direction, width):
    #print("eroding chunk %d,%d in direction %d" % (chunk[0], chunk[1], direction))

    (minX, maxX, minZ, maxZ) = erosionBoundaries(direction, width)
    #print("%s, %d: %d,%d %d,%d" % (Erode.Map[direction], width, minX, maxX, minZ, maxZ))
    
    for x in range(minX, maxX):
        for z in range(minZ, maxZ):
            h = erosionHeight(x, z, direction, width)
            
            if h == MAX_HEIGHT:
                continue
            
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
            if h <= WATER_HEIGHT + 1:
                if h <= WATER_HEIGHT:
                    chunk.Blocks[x, z, h : WATER_HEIGHT + 1] = waterID
                # Turn non-water, non-ice blocks along the shoreline, or under the water, into sand.
                if chunk.Blocks[x, z, h - 1] != waterID \
                        and chunk.Blocks[x, z, h - 1] != iceID:
                    chunk.Blocks[x, z, h - 1] = sandID

def main():
    usage = """
bestofboth --find-edges <path_to_world> --edge-file <edge_file>
bestofboth --smooth <path_to_world> --edge-file <edge_file> [--width <1-16>]
"""
    parser = OptionParser(usage = usage)
    parser.add_option("--find-edges", dest="find_edges",
                    help="world to examine")
    parser.add_option("--smooth", dest="smooth", 
                    help="world to smooth")
    parser.add_option("--edge-file", dest="edge_file", 
                    help="edge file to use")
    parser.add_option("--width", dest="width", 
                    default = "16",
                    help="width of the river")

    (options, args) = parser.parse_args()
    
    if options.find_edges and options.smooth:
        parser.error("--find-edges and --smooth can't be specified " \
            "at the same time. Please run with --find-edges first, " \
            "then use the same --edge-file parameter when running " \
            "--smooth.")
    elif not (options.find_edges or options.smooth):
        parser.error("Must specify --find-edges or --smooth.")
    elif not options.edge_file:
        parser.error("Must specify --edge-file.")
    elif options.find_edges and not os.path.exists(os.path.join(options.find_edges, "level.dat")):
        parser.error("'%s' is not a Minecraft world directory (no " \
            "level.dat file was found)." % (options.find_edges))
    elif options.smooth and not os.path.exists(options.edge_file):
        parser.error("Edge file '%s' does not exist, but the edge " \
            "file must exist when --smooth is specified." \
            % (options.edge_file))
    elif options.width and (int(options.width) < 1 or int(options.width) > 16):
        parser.error("--width must be between 1 and 16 (inclusive)")
        
    
    # Phew! Now that the arguments have been validated...
    if options.find_edges:
        find_edges(options.find_edges, options.edge_file)
    elif options.smooth:
        smooth(options.smooth, options.width, options.edge_file)
    
    return


    if os.path.exists(sys.argv[3]):
        shutil.rmtree(sys.argv[3])
    shutil.copytree(sys.argv[1], sys.argv[3])

    level1 = mclevel.fromFile(sys.argv[3])
    level2 = mclevel.fromFile(sys.argv[2])

    level1bounds = level1.getWorldBounds()
    level2bounds = level2.getWorldBounds()

    print(len(list(level1.allChunks)))
    print(len(list(level2.allChunks)))
    print(level1bounds)
    print(level2bounds)

    print("destination world bounds", level1bounds)
    print("source world bounds", level2bounds)

    """
    # Center world 2 in world 1
    copyTo = [(level1bounds.origin[i] + (level1bounds.size[i] - level2bounds.size[i]) / 2) for i in range(0, 3)]

    print(copyTo)
    shift = [(copyTo[i] - level2bounds.origin[i]) for i in range(0, 3)]
    print("shift: ", shift)
    """

    print("merging worlds...")
    level1.copyBlocksFrom(level2, level2bounds, level2bounds.origin)
    
    for i in range(0, NUM_WORKERS):
        workers.append(Eroder(i, level1, erodeQueue, progressQueue))
        workers[i].start()
    
    numTasks = 0
    level2ChunkSet = set(level2.allChunks)
    for chunk in level2.allChunks:
        numTasks += checkChunk(level2, chunk, erodeQueue)
    
    start = time.time()
    print(start)
    """
    for (chunkCoords, direction) in toErode:
        # Things have been shifted. Make sure we're using new coordinates.
        #chunkCoords = (chunkCoords[0] + shift[X] / 16, chunkCoords[1] + shift[Z] / 16)
        if chunkCoords in level1.allChunks:
            chunk = level1.getChunk(chunkCoords[0], chunkCoords[1])
            erodeQueue.put((chunk, direction))
            numTasks += 1
    """
    
    # One null task is added to the end of the queue for each worker. A
    # worker will stop immediately when it gets a null task from the
    # queue.
    for worker in workers:
        erodeQueue.put((None, None))
    
    for worker in workers:
        worker.join()
    
    while not progressQueue.empty():
        coords = progressQueue.get(True)
        try:
            chunk = level1.getChunk(coords[0], coords[1])
            chunk.chunkChanged()
        except ChunkNotPresent as cnp:
            pass
        
    end = time.time()
    print(end)
    print("elapsed: %f" % (end - start))

    level1.saveInPlace()

if __name__ == "__main__":
    main()
