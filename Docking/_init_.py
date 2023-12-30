from gym.envs.registration import register

register(
    id='fake-v0',
    entry_point='gym_fake.envs:FakeEnv',
)
register(
    id='fake-extrahard-v0',
    entry_point='gym_fake.envs:FakeExtraHardEnv',
)
