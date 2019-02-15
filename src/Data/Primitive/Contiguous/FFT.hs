{-# language BangPatterns        #-}
{-# language LambdaCase          #-}
{-# language NoImplicitPrelude   #-}
{-# language ScopedTypeVariables #-}

-- | This module exposes functions for performing
--   Fast Fourier Transform (FFT) and Inverse Fast Fourier Transform (IFFT)
--   over 'Contiguous' data structures.
module Data.Primitive.Contiguous.FFT
  ( fft
  , ifft
  ) where

import qualified Prelude

import Data.Bool (Bool,otherwise)
import Data.Bits (shiftR,shiftL,(.&.),(.|.))
import Data.Semiring (negate,(+),(*),(-))
import Control.Applicative (pure)
import Control.Monad (when)
import Data.Eq (Eq(..))
import Data.Ord (Ord(..))
import Data.Function (($),(.))
import Data.Complex (Complex(..),conjugate)
import GHC.Exts
import GHC.Real ((/))
import Control.Monad.ST (ST,runST)
import Data.Primitive.Contiguous (Contiguous,Element,Mutable)
import qualified Data.Primitive.Contiguous as Contiguous

{-# RULES 
"fft/ifft" forall x. fft (ifft x) = x
"ifft/fft" forall x. ifft (fft x) = x
  #-}

-- | Radix-2 decimation-in-time fast Fourier Transform.
--   The given array must have a length that is a power of two.
fft :: forall arr. (Contiguous arr, Element arr (Complex Double))
  => arr (Complex Double)
  -> arr (Complex Double)
{-# inlinable [1] fft #-}
fft arr = if arrOK arr
  then runST $ do {
      marr <- copyWhole arr
    ; mfft marr
    ; Contiguous.unsafeFreeze marr
  }
  else Prelude.error "Data.Primitive.Contiguous.FFT.fft: bad array length"

-- | Inverse fast Fourier transform.
ifft :: forall arr. (Contiguous arr, Element arr (Complex Double))
  => arr (Complex Double)
  -> arr (Complex Double)
{-# inlinable [1] ifft #-}
ifft arr = if arrOK arr
  then
    let lenComplex = intToComplexDouble (Contiguous.size arr)
    in cmap ((/lenComplex) . conjugate) . fft . cmap conjugate $ arr
  else Prelude.error "Data.Primitive.Contiguous.FFT.ifft: bad vector length"

copyWhole :: forall arr s a. (Contiguous arr, Element arr a)
  => arr a
  -> ST s (Mutable arr s a)
{-# inline copyWhole #-}
copyWhole arr = do
  let len = Contiguous.size arr
  marr <- Contiguous.new len
  Contiguous.copy marr 0 arr 0 len
  pure marr

arrOK :: forall arr a. (Contiguous arr, Element arr a)
  => arr a
  -> Bool
{-# inline arrOK #-}
arrOK arr =
  let n = Contiguous.size arr
  in (1 `shiftL` log2 n) == n
 
-- array length must be power of two. This property is not checked 
mfft :: forall arr s. (Contiguous arr, Element arr (Complex Double))
  => Mutable arr s (Complex Double)
  -> ST s ()
mfft mut = do {
    len <- Contiguous.sizeMutable mut 
  ; let bitReverse !i !j = do {
          ; if i == len - 1
              then stage 0 1
              else do {
                  when (i < j) $ swap mut i j
                ; let inner k l = if k <= l
                        then inner (k `shiftR` 1) (l - k)
                        else bitReverse (i + 1) (l + k)
                ; inner (len `shiftR` 1) j 
              }
        }
        stage l l1 = if l == (log2 len)
          then pure ()
          else do {
              let !l2 = l1 `shiftL` 1
                  !e = (negate twoPi) / (intToDouble l2)
                  flight j !a = if j == l1
                    then stage (l + 1) l2
                    else do {
                        let butterfly i = if i >= len
                              then flight (j + 1) (a + e)
                              else do {
                                  let i1 = i + l1
                                ; xi1 :+ yi1 <- Contiguous.read mut i1
                                ; let !c = Prelude.cos a
                                      !s = Prelude.sin a
                                      d = (c*xi1 - s*yi1) :+ (s*xi1 + c*yi1)
                                ; ci <- Contiguous.read mut i
                                ; Contiguous.write mut i1 (ci - d)
                                ; Contiguous.write mut i (ci + d)
                                ; butterfly (i + l2)
                              }
                      ; butterfly j
                    }
            ; flight 0 0
         }
  ; bitReverse 0 0
}

-- wildcard cases should never happen. if they do, really bad things will happen.
b,s :: Int -> Int
b = \case { 0 -> 0x02; 1 -> 0x0c; 2 -> 0xf0; 3 -> 0xff00; 4 -> wordToInt 0xffff0000; 5 -> wordToInt 0xffffffff00000000; _ -> 0; }
s = \case { 0 -> 1; 1 -> 2; 2 -> 4; 3 -> 8; 4 -> 16; 5 -> 32; _ -> 0; }
{-# inline b #-}
{-# inline s #-}

log2 :: Int -> Int
log2 v0 = if v0 <= 0
  then Prelude.error $ "Data.Primitive.Contiguous.FFT: nonpositive input, got " Prelude.++ Prelude.show v0
  else go 5 0 v0
  where
    go !i !r !v
      | i == -1 = r
      | v .&. b i /= 0 =
          let si = s i
          in go (i - 1) (r .|. si) (v `shiftR` si)
      | otherwise = go (i - 1) r v

swap :: forall arr s x. (Contiguous arr, Element arr x)
  => Mutable arr s x
  -> Int
  -> Int
  -> ST s ()
{-# inline swap #-}
swap mut i j = do
  atI <- Contiguous.read mut i
  atJ <- Contiguous.read mut j
  Contiguous.write mut i atJ
  Contiguous.write mut j atI

twoPi :: Double
{-# inline twoPi #-}
twoPi = 6.283185307179586

intToDouble :: Int -> Double
{-# inline intToDouble #-}
intToDouble = Prelude.fromIntegral

wordToInt :: Word -> Int
{-# inline wordToInt #-}
wordToInt = Prelude.fromIntegral

intToComplexDouble :: Int -> Complex Double
{-# inline intToComplexDouble #-}
intToComplexDouble = Prelude.fromIntegral

cmap :: (Contiguous arr, Element arr (Complex Double))
  => (Complex Double -> Complex Double)
  -> arr (Complex Double)
  -> arr (Complex Double)
{-# inline cmap #-}
cmap = Contiguous.map
