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
    
    def run(self, level, erosionWidth = 8, waterWidth = 2, maxPreservedHeight = 90):
        raise Exception("not implemented")
    
    def erode(self, chunk, x, z, h):
        #print("setting (%d,%d) to h=%d" % (x, z, h))
        
        if h == MAX_HEIGHT:
            return
        
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
        
class CornerErosionTask(ErosionTask):
    def __init__(self, cornerDirection, cornerPosX, cornerPosZ):
        ErosionTask.__init__(self, cornerPosX, cornerPosZ)
        self.cornerDirection = cornerDirection
        
    def __repr__(self):
        return "corner %-2s %d %d" % (Erode.Map[self.cornerDirection], self.posX, self.posZ)
        
    def run(self, level, erosionWidth = 8, waterWidth = 3, maxPreservedHeight = 75):
        chunksToEdit = {}
        
        # First, make sure all the chunks we need are present.
        try:
            for cx in range(-1, 1):
                for cz in range(-1, 1):
                    chunk = level.getChunk(self.posX / 16 + cx, self.posZ / 16 + cz)
                    chunksToEdit[(cx, cz)] = chunk
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
                print("trying to change CORNER chunk %d,%d.. hp %d,%d" % (self.posX + cx, self.posZ + cz, highPointX, highPointZ))
                for x in range(-8 * cx, -8 * (cx - 1)):
                    for z in range(-8 * cz, -8 * (cz - 1)):
                        dx = x - highPointX
                        dz = z - highPointZ
                        
                        distanceFromCenter = math.sqrt(dx * dx + dz * dz)
                        distanceFromEdge = 8 - distanceFromCenter
        
                        if distanceFromEdge < 0:
                            distanceFromEdge *= -1
                        
                        # We're far enough from the edge that no erosion should occur.
                        if distanceFromEdge > erosionWidth:
                            continue
                        
                        distanceFromRiver = distanceFromEdge - waterWidth
                        if distanceFromRiver < 0:
                            distanceFromRiver = 0
                        
                        # (distanceFromEdge / halfErosionWidth) is on the interval [0..1],
                        # so h will be on the interval [2^0..2^1], or [1..2].
                        h = 2 ** (float(distanceFromRiver) / float(erosionWidth - waterWidth))
                        
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
                        
                        self.erode(chunk, x, z, int(h))
                        chunkChanged = True
                print("changed a corner chunk (%s %d,%d)" % (Erode.Map[self.cornerDirection], self.posX + cx, self.posZ + cz))
                if chunkChanged:
                    chunk.chunkChanged()
        return True
        
class EdgeErosionTask(ErosionTask):
    def __init__(self, edgeDirection, edgePosX, edgePosZ):
        ErosionTask.__init__(self, edgePosX, edgePosZ)
        self.edgeDirection = edgeDirection
        
    def __repr__(self):
        return "edge   %-2s %d %d" % (Erode.Map[self.edgeDirection], self.posX, self.posZ)
        
    def run(self, level, erosionWidth = 8, waterWidth = 3, maxPreservedHeight = 75):
        chunksToEdit = {}
        
        # First, make sure all the chunks we need are present.
        try:
            for cx in range(-1, 1):
                for cz in range(-1, 1):
                    chunk = level.getChunk(self.posX / 16 + cx, self.posZ / 16 + cz)
                    chunksToEdit[(cx, cz)] = chunk
        except ChunkNotPresent:
            return False
        
        for cx in range(-1, 1):
            for cz in range(-1, 1):
                chunkChanged = False
                xMin = -8 * cx
                xMax = -8 * (cx - 1)
                zMin = -8 * cz
                zMax = -8 * (cz - 1)
                print("trying to change chunk %d,%d from x=%d..%d, z=%d..%d" % (self.posX + cx, self.posZ + cz, xMin, xMax, zMin, zMax))
                chunk = chunksToEdit[(cx, cz)]
                for x in range(xMin, xMax):
                    for z in range(zMin, zMax):
                        if self.edgeDirection in [Erode.T, Erode.B]:
                            distanceFromCenter = abs(8 - z)
                        elif self.edgeDirection in [Erode.L, Erode.R]:
                            distanceFromCenter = abs(8 - x)
                        else:
                            raise Exception("unrecognized edge direction %d (%s)" % (self.edgeDirection, Erode.Map[self.edgeDirection]))
                            
                        distanceFromEdge = 8 - distanceFromCenter
        
                        if distanceFromEdge < 0:
                            distanceFromEdge *= -1
                        
                        # We're far enough from the edge that no erosion should occur.
                        if distanceFromEdge > erosionWidth:
                            continue
                        
                        distanceFromRiver = distanceFromEdge - waterWidth
                        if distanceFromRiver < 0:
                            distanceFromRiver = 0
                        
                        distanceFromRiver = distanceFromEdge - waterWidth
                        
                        # (distanceFromEdge / halfErosionWidth) is on the interval [0..1],
                        # so h will be on the interval [2^0..2^1], or [1..2].
                        h = 2 ** (float(distanceFromRiver) / float(erosionWidth - waterWidth))
                        
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
                        
                        self.erode(chunk, x, z, int(h))
                        chunkChanged = True
                if chunkChanged:
                    print("changed an edge chunk (%s %d,%d)" % (Erode.Map[self.edgeDirection], self.posX, self.posZ))
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
    Map = [ "TL", "T", "TR", "L", "R", "BL", "B", "BR",
            "NTL", "NTR", "NBL", "NBR" ]
    Corner = 90
    Edge = 91
    
    TL  = 0
    T   = 1
    TR  = 2
    L   = 3 # toward the left
    R   = 4 # toward the right
    BL  = 5 # toward the bottom-left 
    B   = 6
    BR  = 7

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
                addEdge(coords, toErode, Erode.T)
            
            if not (neighbors[B] or neighbors[BL]):
                coordsBelow = (coords[0], coords[1] + 1)
                addEdge(coordsBelow, toErode, Erode.T)
        
        if neighbors[T]:
            if not (neighbors[L] or neighbors[TL]):
                addEdge(coords, toErode, Erode.L)
        
            if not (neighbors[R] or neighbors[TR]):
                coordsRight = (coords[0] + 1, coords[1])
                addEdge(coordsRight, toErode, Erode.L)
        
    for task in toErode:
        erodeQueue.put(task)
        tasksAdded += 1
    
    return tasksAdded

# erosionWidth: half of the total width of the erosion area, from the
#   middle of the river to one edge
# riverWidth: half of the width, in blocks, of the river running through
#   the erosion area
# maxPreservedHeight: Within the erosion area, everything above this
#   height will be eroded away, even on the edges of the erosion area.
#   There will be a slope between the edges and the river.

def erosionHeight(x, z, erodeType, xCenter, zCenter, erosionWidth = 8, riverWidth = 4, maxPreservedHeight = 90):
    
    h = MAX_HEIGHT
    if erodeType == Erode.Corner:
        dx = xCenter - x - 0.5
        dz = zCenter - z - 0.5
        distanceFromCenter = math.sqrt(dx * dx + dz * dz)
        distanceFromEdge = 8 - distanceFromCenter
        
        if distanceFromEdge < 0:
            distanceFromEdge *= -1
        
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
    
    return int(h)
            
def erosionHeightOld(x, z, direction, erosionWidth = 8, riverWidth = 1, maxPreservedHeight = 90):
    
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

def erode(chunk, erodeType, xCenter, zCenter, xMin, xMax, zMin, zMax, width):
    #print("eroding x=%d..%d, z=%d..%d" % (xMin, xMax, zMin, zMax))
    for x in range(xMin, xMax):
        for z in range(zMin, zMax):
            h = erosionHeight(x, z, erodeType, xCenter, zCenter, width)
            
            #if (xCenter == 8 and zCenter == 8):
            #    print("eroding w%d, %d,%d... h=%d" % (width, x, z, h))
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

def erode_old(chunk, direction, width):
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
    parser.add_option("--width", dest="width", 
                    default = "16",
                    help="width of the river")

    (options, args) = parser.parse_args()
    
    worldDir = options.find_edges or options.smooth
    
    if options.find_edges and options.smooth:
        parser.error("--find-edges and --smooth can't be specified " \
            "at the same time. Please run with --find-edges first, " \
            "then use the same --edge-file parameter when running " \
            "--smooth.")
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
