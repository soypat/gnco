package main

import (
	"context"
	"math/rand"

	"github.com/soypat/mu8"
	"github.com/soypat/mu8/genes"
	"github.com/soypat/mu8/genetic"
)

func main() {

}

func run() error {

	genetic.NewIslands(8, make([]*GeneticRocket, 120), rand.NewSource(1), newIndividual)
}

func newIndividual() *GeneticRocket {
	return &GeneticRocket{}
}

var _ mu8.Genome = (*GeneticRocket)(nil)

type GeneticRocket struct {
	genes.ConstrainedFloatGrad
}

func (gr *GeneticRocket) Simulate(ctx context.Context) (fitness float64) {
	return
}

func (gr *GeneticRocket) GetGene(i int) mu8.Gene {
	return nil
}

func (gr *GeneticRocket) Len() int {
	return 0
}
