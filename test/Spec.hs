{-# language TemplateHaskell #-}

{-# options_ghc -O0 #-}

module Main (main) where

import Data.Complex (Complex(..))

import Data.Primitive.PrimArray (PrimArray)
import qualified Data.Primitive.PrimArray as PA

import Hedgehog
import qualified Hedgehog.Gen as Gen
import qualified Hedgehog.Range as Range

import Data.Primitive.Instances ()

import qualified Data.Primitive.Contiguous.FFT as CF
import qualified Numeric.MathFunctions.Comparison as MF

epsilon :: Double
epsilon = 0.0000000003

main :: IO Bool
main = checkParallel $$(discover)

prop_close :: Property
prop_close = property $ do
  x <- forAll genPrimArray
  let fft' = fft x
      {-# noinline fft' #-}
  let ifft' = ifft fft'
      {-# noinline ifft' #-}
  pwithin ifft' x === True

fft,ifft :: PrimArray (Complex Double) -> PrimArray (Complex Double)
fft = CF.fft
ifft = CF.ifft

pwithin :: PrimArray (Complex Double) -> PrimArray (Complex Double) -> Bool
pwithin xs ys = all (== True)
  $ zipWith within (PA.primArrayToList xs) (PA.primArrayToList ys)

within :: Complex Double -> Complex Double -> Bool
within (x :+ y) (x' :+ y') = abs (x - x') <= epsilon
  && abs (y - y') <= epsilon

genDouble :: Gen Double
genDouble = Gen.double (Range.linearFrac 5 200)

genComplex :: Gen (Complex Double)
genComplex = (:+) <$> genDouble <*> genDouble

genPrimArray :: Gen (PrimArray (Complex Double))
genPrimArray = PA.primArrayFromList <$> Gen.list (Range.singleton 8) genComplex
