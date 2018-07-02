{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE NoImplicitPrelude   #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Data.Primitive.Contiguous.FFT
  ( dft
  , idft
  , dftMutable
  , subDFT
  ) where

import qualified Prelude

import Data.Eq (Eq((==)))
import Data.Function (($))
import Control.Monad
import Data.Ord
import Control.Monad.ST
import Data.Complex hiding (cis)
import qualified Data.Complex as C
import Data.Primitive.Contiguous
import GHC.Num (Num(..))
import GHC.Float
import GHC.Real
import GHC.Exts (Int)

cis :: Floating a => a -> a -> Complex a
cis k n = C.cis (2 * pi * k / n)
{-# INLINE cis #-}

mkComplex :: x -> x -> Complex x
mkComplex !r !i = r :+ i
{-# INLINE mkComplex #-}

clone :: Contiguous arr => Element arr b => Mutable arr s b -> ST s (Mutable arr s b)
clone !m = sizeMutable m >>= \x -> cloneMutable m 0 (x - 1)

dftMutable :: forall arr x s. (RealFloat x, Contiguous arr, Element arr (Complex x))
     => Mutable arr s (Complex x)
     -> ST s (Mutable arr s (Complex x))
dftMutable !mut = do
  !sz <- sizeMutable mut
      
  let getII !ix = (ix + sz `Prelude.div` 2) `Prelude.mod` sz 
      go :: Int -- ^ i value
         -> Int -- ^ j value
         -> Complex x -- ^ accumulator
         -> ST s ()
      go !i !j !acc = if i == sz then return () else if j < sz
        then do
          let !jj = getII j
          atJJ@(r :+ _) <- read mut jj
          let real, imag, same :: x
              !same = (-2) * pi * (fromIntegral (i * j)) / (fromIntegral sz)
              !real = r * cos same
              !imag = r * sin same
              !val  = acc + mkComplex real imag
          go i (j + 1) val
        else do
          let !ii = getII i
          !_ <- write mut ii acc :: ST s ()
          go (i + 1) 0 0

  !_ <- go 0 0 0

  return mut

dft :: forall arr x. (RealFloat x, Contiguous arr, Element arr x, Element arr (Complex x))
     => arr x
     -> arr (Complex x)
dft !a = runST $ dftInternal a

-- | not in-place, also very inefficient. currently /O(n^2)/
dftInternal :: forall arr x s. (RealFloat x, Contiguous arr, Element arr x, Element arr (Complex x))
  => arr x
  -> ST s (arr (Complex x))
dftInternal !a = do
  let !sz = size a
      getII !ix = (ix + sz `Prelude.div` 2) `Prelude.mod` sz 
  
  !mut <- new sz :: ST s (Mutable arr s (Complex x))
 
  let go :: Int -- ^ i value
         -> Int -- ^ j value
         -> Complex x -- ^ accumulator
         -> ST s ()
      go !i !j !acc = if i == sz then return () else if j < sz
        then do
          let !jj = getII j
              !atJJ = index a jj
              real, imag, same :: x
              !same = (-2) * pi * (fromIntegral (i * j)) / (fromIntegral sz)
              !real = atJJ * cos same
              !imag = atJJ * sin same
              !val  = acc + mkComplex real imag
          go i (j + 1) val
        else do
          let !ii = getII i
          !_ <- write mut ii acc :: ST s ()
          go (i + 1) 0 0

  !_ <- go 0 0 0

  unsafeFreeze mut

idft :: forall arr x. (RealFloat x, Contiguous arr, Element arr x, Element arr (Complex x))
  => arr (Complex x)
  -> arr x
idft !a = runST $ idftInternal a

-- | not in-place, also very inefficient. currently /O(n^2)/
idftInternal :: forall arr x s. (RealFloat x, Contiguous arr, Element arr x, Element arr (Complex x))
  => arr (Complex x)
  -> ST s (arr x)
idftInternal !a = do
  let !sz = size a
      getII !ix = (ix + sz `Prelude.div` 2) `Prelude.mod` sz

  !mut <- new sz :: ST s (Mutable arr s x)
  
  let go :: Int
         -> Int
         -> x
         -> ST s ()
      go !i !j !acc = if i == sz then return () else if j < sz
        then do
          let !jj = getII j
              !atJJ@(real :+ imag) = index a jj
              !sCount = fromIntegral sz
              !same = (-2) * pi * (fromIntegral (i * j)) / sCount
              !val = (real * cos same + imag * sin same) / sCount
          go i (j + 1) val
        else do
          let !ii = getII i 
          !_ <- write mut ii acc :: ST s ()
          go (i + 1) 0 0

  !_ <- go 0 0 0

  unsafeFreeze mut

-- | FIXME: doc
--   FIXME: name
--   FIXME: source
--
--   Given a signal size, previous window, transform of previous window, and the newest value,
--   compute the transform of the new window (which is just a shifted version of the previous window)
--   in /O(n)/ time, in-place
subDFT :: forall arr x s. (RealFloat x, Contiguous arr, Element arr x, Element arr (Complex x))
  => Int   -- ^ N, signal size
  -> Mutable arr s (Complex x) -- ^ x1, original window
  -> Complex x -- ^ newest complex value
  -> Mutable arr s (Complex x) -- ^ f1, previous transform
  -> ST s (Mutable arr s (Complex x)) -- ^ f2, new transform
subDFT n x1 x2_N_1 f1 = do
  let !sz = fromIntegral n :: x

  !l <- sizeMutable f1
  !x1_0 <- read x1 0 :: ST s (Complex x)
 
  let go :: Int -> ST s ()
      go !ix = if ix < l
        then do
          f1_k <- read f1 ix
          let foo' = cis (fromIntegral ix) sz
              res  = f1_k + x2_N_1 + x1_0
              fin  = foo' * res
          !_ <- write f1 ix fin
          go (ix + 1)
        else return ()
  go 0
  return f1
  
