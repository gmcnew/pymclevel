import math
import Queue
import shutil
import sys
import threading

import materials
import mclevel

X = 0
Y = 1
Z = 2

CORNER_TOP_LEFT = 0
CORNER_TOP_RIGHT = 1
CORNER_BOTTOM_LEFT = 2
CORNER_BOTTOM_RIGHT = 3

EDGE_TOP = 0
EDGE_BOTTOM = 1
EDGE_LEFT = 2
EDGE_RIGHT = 3

class Erode:
    TL = 0
    T  = 1
    TR = 2
    L  = 3
    R  = 4
    BL = 5
    B  = 6
    BR = 7
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

leavesID = materials.materials.Leaves.ID
pineLeavesID = materials.materials.PineLeaves.ID
birchLeavesID = materials.materials.BirchLeaves.ID

class Eroder(threading.Thread):
    def __init__(self, ID, erodeQueue, progressQueue, stopSignal):
        threading.Thread.__init__(self)
        self.ID = ID
        self.erodeQueue = erodeQueue
        self.progressQueue = progressQueue
        self.stopSignal = stopSignal
    def run(self):
        while True:
            # Wait forever for a new task. Null tasks (where chunk is
            # None) will be added to the queue when there are no more
            # work items, so this is guaranteed to return an element
            # eventually.
            (chunk, direction) = self.erodeQueue.get(True)
            if chunk is None:
                print("worker %d is returning" % self.ID)
                break
            erode(chunk, direction)
            self.progressQueue.put(chunk)
        

def checkEdgesAndCorners(levelChunkSet, chunk, toErode):
    aroundMe = [(-1, -1), (0, -1), (1, -1),
                (-1,  0),          (1,  0),
                (-1,  1), (0,  1), (1,  1)]
    
    neighbors = [True] * 8
    
    onPerimeter = False
    
    for i in range(len(aroundMe)):
        a = aroundMe[i]
        if (chunk[0] + a[0], chunk[1] + a[1]) not in levelChunkSet:
            onPerimeter = True
            neighbors[i] = False
    
    if not onPerimeter:
        return
    
    # Top-left corner
    if (not neighbors[TL]) and (not neighbors[T]) and (not neighbors[L]):
        toErode.append((chunk, Erode.TL))
        toErode.append(((chunk[0] - 1, chunk[1] - 1), Erode.NBR))
    
    # Top-right corner
    if (not neighbors[T]) and (not neighbors[TR]) and (not neighbors[R]):
        toErode.append((chunk, Erode.TR))
        toErode.append(((chunk[0] + 1, chunk[1] - 1), Erode.NBL))
    
    # Bottom-right corner
    if (not neighbors[R]) and (not neighbors[BR]) and (not neighbors[B]):
        toErode.append((chunk, Erode.BR))
        toErode.append(((chunk[0] + 1, chunk[1] + 1), Erode.NTL))
    
    # Bottom-left corner
    if (not neighbors[B]) and (not neighbors[BL]) and (not neighbors[L]):
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
    
    if neighbors[T] and neighbors[L] and not neighbors[TL]:
        toErode.append((chunk, Erode.NTL))
        toErode.append(((chunk[0] - 1, chunk[1] - 1), Erode.BR))
    
    if neighbors[T] and neighbors[R] and not neighbors[TR]:
        toErode.append((chunk, Erode.NTR))
        toErode.append(((chunk[0] + 1, chunk[1] - 1), Erode.BL))
    
    if neighbors[B] and neighbors[R] and not neighbors[BR]:
        toErode.append((chunk, Erode.NBR))
        toErode.append(((chunk[0] + 1, chunk[1] + 1), Erode.TL))
    
    if neighbors[B] and neighbors[L] and not neighbors[BL]:
        toErode.append((chunk, Erode.NBL))
        toErode.append(((chunk[0] - 1, chunk[1] + 1), Erode.TR))
        
    """
        corners.add((chunk, CORNER_TOP_LEFT))
    if not (neighbors[T] and neighbors[TR] and neighbors[R]):
        corners.add((chunk, CORNER_TOP_RIGHT))
    if not (neighbors[L] and neighbors[BL] and neighbors[B]):
        corners.add((chunk, CORNER_BOTTOM_LEFT))
    if not (neighbors[R] and neighbors[B] and neighbors[BR]):
        corners.add((chunk, CORNER_BOTTOM_RIGHT))
    
    if neighbors[L] and neighbors[R] and \
            (not neighbors[TL]) and (not neighbors[T]) and (not neighbors[TR]):
        edges.add((chunk, EDGE_TOP))
    elif neighbors[T] and neighbors[B] and \
            (not neighbors[TL]) and (not neighbors[L]) and (not neighbors[BL]):
        edges.add((chunk, EDGE_LEFT))
    elif neighbors[T] and neighbors[B] and\
            (not neighbors[TR]) and (not neighbors[R]) and (not neighbors[BR]):
        edges.add((chunk, EDGE_RIGHT))
    if not (neighbors[4]):
        edges.add((chunk, EDGE_RIGHT))
    if not (neighbors[6]):
        edges.add((chunk, EDGE_BOTTOM))
    """

def isPerimeterChunk(level, chunk):
    aroundMe = [(-1, -1), (0, -1), (1, -1),
                (-1,  0),          (1,  0),
                (-1,  1), (0,  1), (1,  1)]
    for a in aroundMe:
        if (chunk[0] + a[0], chunk[1] + a[1]) not in level.allChunks:
            return True
    return False

WATER_HEIGHT = 63

def getErodedHeight(x, z, direction):
    lowPoint = None
    highPoint = None
    isCorner = True
    if direction == Erode.TL:
        highPoint = (8, 8)
    elif direction == Erode.TR:
        highPoint = (8, 8)
    elif direction == Erode.BL:
        highPoint = (8, 8)
    elif direction == Erode.BR:
        highPoint = (8, 8)
    elif direction == Erode.NTL:
        lowPoint = (0, 0)
    elif direction == Erode.NTR:
        lowPoint = (16, 0)
    elif direction == Erode.NBL:
        lowPoint = (0, 16)
    elif direction == Erode.NBR:
        lowPoint = (16, 16)
    else:
        isCorner = False
        
    
    if highPoint:
        dx = highPoint[0] - x - 0.5
        dz = highPoint[1] - z - 0.5
        distanceFromCenter = math.sqrt((dx * dx) + (dz * dz))
        distanceFromEdge = 8 - distanceFromCenter
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
        
    h = int(WATER_HEIGHT + (2 ** distanceFromEdge) - 1)
    
    return h

def erosionBoundaries(direction):
    if direction == Erode.TL or direction == Erode.NTL:
        return (0, 8, 0, 8)
    elif direction == Erode.TR or direction == Erode.NTR:
        return (8, 16, 0, 8)
    elif direction == Erode.BR or direction == Erode.NBR:
        return (8, 16, 8, 16)
    elif direction == Erode.BL or direction == Erode.NBL:
        return (0, 8, 8, 16)
    elif direction == Erode.T:
        return (0, 16, 0, 8)
    elif direction == Erode.B:
        return (0, 16, 8, 16)
    elif direction == Erode.L:
        return (0, 8, 0, 16)
    elif direction == Erode.R:
        return (8, 16, 0, 16)

def erode(chunk, direction, distance = 8):
    #print("eroding chunk %d,%d in direction %d" % (chunk[0], chunk[1], direction))

    (minX, maxX, minZ, maxZ) = erosionBoundaries(direction)
    
    for x in range(minX, maxX):
        for z in range(minZ, maxZ):
            h = getErodedHeight(x, z, direction)
            
            # Turn everything in this vertical column into air, but
            # leave leaves alone to avoid weird-looking half-trees.
            for h2 in range(h, 128):
                if chunk.Blocks[x, z, h2] not in [leavesID, pineLeavesID, birchLeavesID]:
                    chunk.Blocks[x, z, h2] = airID
            if h <= WATER_HEIGHT + 1:
                if h == WATER_HEIGHT:
                    chunk.Blocks[x, z, h] = waterID
                # Turn non-water, non-ice blocks along the shoreline, or under the water, into sand.
                if chunk.Blocks[x, z, h - 1] != waterID \
                        and chunk.Blocks[x, z, h - 1] != iceID:
                    chunk.Blocks[x, z, h - 1] = sandID

NUM_WORKERS = 2

def main():
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

    # Center world 2 in world 1

    #   l2o + (l1c - l2c)
    # = l2o + ((l1o + l1s/2) - (l2o + l2s/2))
    # = l2o + l1o + l1s/2 - l2o - l2s/2
    # = l1o + (l1s - l2s) / 2
    """
    level1center = [(level1bounds.origin[i] + (level1bounds.size[i] / 2)) for i in range(0, 3)]
    level1center[Y] = 0

    level2center = [(level2bounds.origin[i] + (level2bounds.size[i] / 2)) for i in range(0, 3)]
    level2center[Y] = 0

    print("level 1 center", level1center)
    print("level 2 center", level2center)
    """

    #copyTo = [(level2bounds.origin[i] + (level1center[i] - level2center[i])) for i in range(0, 3)]
    copyTo = [(level1bounds.origin[i] + (level1bounds.size[i] - level2bounds.size[i]) / 2) for i in range(0, 3)]
    #copyTo[1] = 0
    #copyTo = [-240, 0, -128]#level1bounds.origin

    print(copyTo)
    shift = [(copyTo[i] - level2bounds.origin[i]) for i in range(0, 3)]
    print("shift: ", shift)

    print("merging worlds...")
    level1.copyBlocksFrom(level2, level2bounds, copyTo)

    edges = set()
    corners = set()
    toErode = []
    
    level2ChunkSet = set(level2.allChunks)
    for chunk in level2.allChunks:
        checkEdgesAndCorners(level2ChunkSet, chunk, toErode)
        #if isPerimeterChunk(level2, chunk):
        #    perimeterChunks.add(chunk)
        
    workers = []
    erodeQueue = Queue.Queue()
    progressQueue = Queue.Queue()
    stopSignal = threading.Event()
    for i in range(0, NUM_WORKERS):
        workers.append(Eroder(i, erodeQueue, progressQueue, stopSignal))
        workers[i].start()
    
    chunksChanged = set()
    import time
    start = time.time()
    print(start)
    eroded = 0
    lastProgress = 0
    
    numTasks = 0
    
    for (chunkCoords, direction) in toErode:
        # Things have been shifted. Make sure we're using new coordinates.
        chunkCoords = (chunkCoords[0] + shift[X] / 16, chunkCoords[1] + shift[Z] / 16)
        if chunkCoords in level1.allChunks:
            chunk = level1.getChunk(chunkCoords[0], chunkCoords[1])
            erodeQueue.put((chunk, direction))
            numTasks += 1

            """
            toErode.put(
            erode(level1, chunkCoords, direction)
            eroded += 1
            progress = int(eroded * 100 / len(toErode))
            if progress != lastProgress:
                lastProgress = progress
                print("progress: %d%%" % (progress))
            chunksChanged.add(chunkCoords)
            """
    
    # One null task is added to the end of the queue for each worker. A
    # worker will stop immediately when it gets a null task from the
    # queue.
    for worker in workers:
        erodeQueue.put((None, None))
    
    tasksFinished = 0
    while tasksFinished < numTasks:
        processedChunk = progressQueue.get(True)
        processedChunk.chunkChanged()
        tasksFinished += 1
        print("progress: %d" % (tasksFinished * 100 / numTasks))
        
    end = time.time()
    print(end)
    print("elapsed: %f" % (end - start))
    
    """
    modified = 0
    for c in chunksChanged:
        chunk = level1.getChunk(c[0], c[1])
        chunk.chunkChanged()
        modified += 1
    print("eroded %d chunks" % modified)
    """

    level1.saveInPlace()

if __name__ == "__main__":
    main()
